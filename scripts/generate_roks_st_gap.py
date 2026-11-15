"""Generate immutable frozen-ROKS probes from preserved singlet/triplet donors.

Run as a development module or as a standalone script beside
roks_st_gap_common.py. No calculations are submitted by this command.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
try:
    from . import roks_st_gap_common as c
except ImportError:
    import roks_st_gap_common as c

MOLECULES = ['Heptazine', 'Cyclazine', *[f'Molecule{i}' for i in range(3, 9)]]


def section(text, name):
    match = re.search(rf'(?ms)^([ \t]*)&{name}\b.*?^\1&END {name}[ \t]*$', text)
    if not match:
        raise ValueError(f'Missing {name}')
    return match[0]


def prepare(source, destination, pairs, molecules):
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise ValueError('Destination must be new; existing snapshots are immutable')
    donors = sorted({d for e, d in pairs})
    manifest = dict(protocol='first-pre-update-explicit-determinant-v1',
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    basis='aug-cc-pVTZ', molecules=molecules, pairs=pairs,
                    qualification='Frozen diagnostic, not an optimized UKS state or a spectroscopic assignment.',
                    files={}, trials=[])
    for molecule in molecules:
        for donor in donors:
            src = source / donor / molecule
            dst = destination / 'sources' / donor / molecule
            dst.mkdir(parents=True)
            states = {}
            for stage in ('singlet', 'triplet'):
                matches = list(src.glob(f'*_{stage.upper()}-RESTART.wfn'))
                if len(matches) != 1:
                    raise ValueError(f'Ambiguous/missing donor checkpoint: {src}/{stage}')
                for path in [matches[0], src / f'{stage}.inp', src / f'{stage}.out']:
                    shutil.copy2(path, dst / path.name)
                states[stage] = c.read_wfn(matches[0])
            st = (src / 'singlet.inp').read_text()
            tt = (src / 'triplet.inp').read_text()
            c.check_basis(st)
            c.check_basis(tt)
            c.check_poisson_pair(st, tt)
            for name in ('SUBSYS', 'QS', 'MGRID', 'XC'):
                if section(st, name) != section(tt, name):
                    raise ValueError(f'{src}: S/T differ in {name}')
            metric, metric_info = c.recover_overlap(states['triplet'])
            orth = c.orthonormality_error(states['singlet'], metric)
            for evaluator, orbital in pairs:
                if orbital != donor:
                    continue
                template = source / evaluator / molecule / 'singlet.inp'
                text = template.read_text()
                for name in ('SUBSYS', 'MGRID', 'QS'):
                    if section(text, name) != section(st, name):
                        raise ValueError(f'{molecule}: incompatible {name} for {evaluator}@{donor}')
                for stage, config in [('singlet', 'M'), ('singlet', 'T'), ('triplet', 'T')]:
                    trial = f'{molecule}__{donor}_{stage}__{evaluator}_{config}'
                    directory = destination / 'trials' / trial
                    directory.mkdir(parents=True)
                    det = c.explicit_determinant(states[stage], config)
                    (directory / 'seed.wfn').write_bytes(c.wfn_bytes(det))
                    (directory / 'probe.inp').write_text(c.probe_input(text, config))
                    wfn = next(dst.glob(f'*_{stage.upper()}-RESTART.wfn'))
                    row = dict(id=trial, molecule=molecule, donor=donor,
                               evaluator=evaluator, source_stage=stage, configuration=config,
                               source_wfn=str(wfn.relative_to(destination)),
                               source_wfn_sha256=c.sha(wfn),
                               seed_sha256=c.sha(directory / 'seed.wfn'),
                               input_sha256=c.sha(directory / 'probe.inp'),
                               template_sha256=c.sha(template),
                               source_orthonormality_error=orth,
                               seed_orthonormality_error=c.orthonormality_error(det, metric),
                               independent_triplet_metric=metric_info)
                    (directory / 'trial.json').write_text(json.dumps(row, indent=2)+'\n')
                    manifest['trials'].append(row)
    # A cross-evaluator need not also be an orbital donor.
    for evaluator in {e for e, d in pairs}:
        for molecule in molecules:
            target = destination / 'sources' / evaluator / molecule / 'singlet.inp'
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / evaluator / molecule / 'singlet.inp', target)
    for path in sorted((destination / 'sources').rglob('*')):
        if path.is_file():
            manifest['files'][str(path.relative_to(destination))] = c.sha(path)
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (destination / 'trials.txt').write_text(''.join(r['id']+'\n' for r in manifest['trials']))
    print(f"Prepared {len(manifest['trials'])} probes at {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Root with functional/molecule donor directories')
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--pairs', nargs='+', default=['r2scan@r2scan','pbe0@pbe0'], help='Evaluator@orbital-donor pairs')
    parser.add_argument('--molecules', nargs='+', choices=MOLECULES, default=MOLECULES)
    args = parser.parse_args()
    pairs = [p.lower().split('@') for p in args.pairs]
    if any(len(p)!=2 or any(f not in ('pbe','r2scan','pbe0','skala') for f in p) for p in pairs) or len({tuple(p) for p in pairs})!=len(pairs):
        parser.error('Pairs must be unique evaluator@donor combinations')
    prepare(args.source, args.destination, pairs, args.molecules)

if __name__ == '__main__':
    main()
