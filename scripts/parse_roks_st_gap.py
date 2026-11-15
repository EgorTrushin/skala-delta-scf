"""Parse archived frozen ROKS evaluations; standalone beside roks_st_gap_common.py.

Default parsing verifies the stored density-validation records and file hashes.
--validate-densities repeats the numerical check and requires archived WFN files.
No production energies are substituted for frozen components.
"""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re
try:
    from . import roks_st_gap_common as c
except ImportError:
    import roks_st_gap_common as c

ORDER = ['Heptazine', 'Cyclazine', *[f'Molecule{i}' for i in range(3, 9)]]
ENERGY = r'ENERGY\|\s+Total FORCE_EVAL\s+\(\s*QS\s*\) energy \[hartree\]\s+(\S+)'
REPRODUCTION_TOL_EV = 1e-4


def donor_energy(path):
    text = path.read_text()
    if text.count('PROGRAM STARTED AT') != 1 or 'PROGRAM ENDED AT' not in text:
        raise ValueError(f'{path}: incomplete or concatenated donor output')
    if 'SCF run NOT converged' in text or 'SCF run converged in' not in text:
        raise ValueError(f'{path}: unconverged donor')
    if not re.search(r'DFT\|.*\bROKS\s*$', text, re.M):
        raise ValueError(f'{path}: donor is not ROKS')
    values = re.findall(ENERGY, text)
    if len(values)!=1 or not c.np.isfinite(float(values[0])):
        raise ValueError(f'{path}: ambiguous energy')
    return float(values[0])


def parse_campaign(root, validate_densities=False):
    manifest = json.loads((root / 'manifest.json').read_text())
    if manifest['protocol'] != 'first-pre-update-explicit-determinant-v1':
        raise ValueError(f'{root}: unsupported protocol')
    for name, digest in manifest.get('files', {}).items():
        path = root / name
        if path.suffix=='.wfn' and not path.exists() and not validate_densities:
            continue
        if c.sha(path)!=digest:
            raise ValueError(f'{path}: changed source')
    groups = defaultdict(dict)
    for trial in manifest['trials']:
        directory = root / 'trials' / trial['id']
        if 'runtime_record' in trial and c.sha(root / trial['runtime_record'])!=trial['runtime_sha256']:
            raise ValueError(f'{directory}: changed runtime record')
        if validate_densities:
            c.validate(directory)
        v = json.loads((directory / 'validation.json').read_text())
        for key in ('id','molecule','donor','evaluator','source_stage','configuration','input_sha256','seed_sha256','source_wfn_sha256'):
            if v[key] != trial[key]:
                raise ValueError(f'{directory}: {key} differs from manifest')
        density_errors = [*v['density_max_difference'].values(), *v['density_relative_error'].values()]
        orbital_errors = [v['source_orthonormality_error'],v['seed_orthonormality_error']]
        if (not v.get('validated') or not c.np.isfinite(density_errors+orbital_errors).all()
                or max(density_errors)>c.DENSITY_TOL or max(orbital_errors)>c.ORTHO_TOL):
            raise ValueError(f'{directory}: failed density validation')
        for filename, key in [('probe.inp','input_sha256'),('probe.out','output_sha256')]:
            if c.sha(directory / filename)!=v[key]:
                raise ValueError(f'{directory}: changed {filename}')
        if (directory / 'seed.wfn').exists() and c.sha(directory / 'seed.wfn')!=v['seed_sha256']:
            raise ValueError(f'{directory}: changed seed')
        density = directory / 'density.Log.gz'
        if not density.exists():
            density = directory / 'density.Log'
        if density.exists() and c.hashlib.sha256(c.density_bytes(density)).hexdigest()!=v['density_sha256']:
            raise ValueError(f'{directory}: changed density')
        source_dir = root / 'sources' / trial['donor'] / trial['molecule']
        source_wfns = list(source_dir.glob(f"*_{trial['source_stage'].upper()}-RESTART.wfn"))
        if source_wfns and (len(source_wfns)!=1 or c.sha(source_wfns[0])!=trial['source_wfn_sha256']):
            raise ValueError(f'{directory}: donor WFN differs')
        template = root / 'sources' / trial['evaluator'] / trial['molecule'] / 'singlet.inp'
        c.check_basis(template.read_text())
        if c.sha(template)!=trial['template_sha256']:
            raise ValueError(f'{directory}: changed evaluator template')
        if c.probe_input(template.read_text(),trial['configuration']) != (directory / 'probe.inp').read_text():
            raise ValueError(f'{directory}: scientific input differs from documented conversion')
        text = (directory / 'probe.out').read_text()
        steps = re.findall(r'^\s*(\d+)\s+OT\s+\S+\s+\S+\s+\S+\s+\S+\s+([-+\d.E]+)\s+[-+\d.E]+\s*$',text,re.M)
        if len(steps)!=1 or steps[0][0]!='1' or abs(float(steps[0][1])-v['energy_Ha'])>1e-10:
            raise ValueError(f'{directory}: not one validated pre-update energy')
        components = c.energy_components(text, v['energy_Ha'])
        if trial['evaluator']=='pbe0' and components['exact_exchange']==0:
            raise ValueError(f'{directory}: missing hybrid exchange')
        key = (trial['molecule'],trial['donor'],trial['evaluator'])
        config = (trial['source_stage'],trial['configuration'])
        if config in groups[key]:
            raise ValueError(f'{directory}: duplicate probe')
        groups[key][config] = components
    pairs = manifest.get('pairs', [['pbe','pbe'],['skala','pbe'],['pbe','skala'],['skala','skala']])
    expected = {(m,d,e) for m in manifest['molecules'] for e,d in pairs}
    if set(groups)!=expected:
        raise ValueError(f'{root}: incomplete comparison matrix')
    rows, checks = [], []
    for (mol,donor,evaluator), group in groups.items():
        if set(group)!={('singlet','M'),('singlet','T'),('triplet','T')}:
            raise ValueError(f'{mol}: incomplete components')
        m,t,tt = [group[k] for k in [('singlet','M'),('singlet','T'),('triplet','T')]]
        em,ets,ett = m['total'],t['total'],tt['total']
        singlet = 2*em-ets
        common = 2*(em-ets)*c.HARTREE_TO_EV
        correction = (ets-ett)*c.HARTREE_TO_EV
        gap = (singlet-ett)*c.HARTREE_TO_EV
        if abs(gap-common-correction)>1e-10:
            raise ValueError('Gap identity failed')
        xc = 2*sum(m[k]-t[k] for k in ('xc_grid','xc_atomic','exact_exchange'))*c.HARTREE_TO_EV
        check = dict(molecule=mol,orbital_functional=donor,evaluating_functional=evaluator,
                     campaign=root.name,non_xc_residual_eV=common-xc,
                     component_splittings_eV={k:2*(m[k]-t[k])*c.HARTREE_TO_EV for k in m if k!='total'},
                     negative_triplet_correction=correction < -1e-6)
        if abs(common-xc)>1e-5:
            raise ValueError(f'{mol}/{evaluator}@{donor}: non-XC cancellation failed: {common-xc}')
        if donor==evaluator:
            src = root / 'sources' / donor / mol
            es0,et0 = [donor_energy(src / (s+'.out')) for s in ('singlet','triplet')]
            check.update(singlet_reproduction_eV=(singlet-es0)*c.HARTREE_TO_EV,
                         triplet_reproduction_eV=(ett-et0)*c.HARTREE_TO_EV,
                         gap_reproduction_eV=gap-(es0-et0)*c.HARTREE_TO_EV)
            if max(abs(check[k]) for k in ('singlet_reproduction_eV','triplet_reproduction_eV','gap_reproduction_eV')) > REPRODUCTION_TOL_EV:
                raise ValueError(f'Donor reproduction failed: {check}')
        rows.append(dict(molecule=mol,basis='aug-cc-pVTZ',orbital_functional=donor,
                         evaluating_functional=evaluator,mixed_at_singlet_Ha=em,
                         triplet_at_singlet_Ha=ets,triplet_at_triplet_Ha=ett,
                         singlet_Ha=singlet,common_gap_eV=common,
                         triplet_correction_eV=correction,gap_eV=gap,source_campaign=root.name))
        checks.append(check)
    return rows,checks


def export(campaigns, output, validate_densities=False):
    rows,checks=[],[]
    for root in campaigns:
        rr,cc=parse_campaign(root,validate_densities)
        rows.extend(rr);checks.extend(cc)
    keys=[(r['molecule'],r['orbital_functional'],r['evaluating_functional']) for r in rows]
    if len(set(keys))!=len(keys):
        raise ValueError('Duplicate adoption across campaigns')
    rows.sort(key=lambda r:(ORDER.index(r['molecule']),r['orbital_functional'],r['evaluating_functional']))
    summary=[]
    for donor,evaluator in sorted({(r['orbital_functional'],r['evaluating_functional']) for r in rows}):
        rr=[r for r in rows if r['orbital_functional']==donor and r['evaluating_functional']==evaluator]
        summary.append(dict(orbital_functional=donor,evaluating_functional=evaluator,n=len(rr),
                            mean_common_gap_eV=sum(r['common_gap_eV'] for r in rr)/len(rr),
                            mean_gap_eV=sum(r['gap_eV'] for r in rr)/len(rr),
                            common_inversions=sum(r['common_gap_eV']<0 for r in rr),
                            gap_inversions=sum(r['gap_eV']<0 for r in rr),
                            sign_changes=[r['molecule'] for r in rr if (r['common_gap_eV']<0)!=(r['gap_eV']<0)]))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    output.with_suffix('.validation.json').write_text(json.dumps(dict(hartree_to_eV=c.HARTREE_TO_EV,
        reproduction_tolerance_eV=REPRODUCTION_TOL_EV,checks=checks,summary=summary),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--validate-densities',action='store_true')
    a=p.parse_args();export(a.campaign,a.output,a.validate_densities)

if __name__=='__main__':
    main()
