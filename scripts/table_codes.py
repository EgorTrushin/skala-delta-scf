"""Print and export CP2K/PySCF MADs and MaxADs from supplied CSV results."""

import argparse
import hashlib
import math
from pathlib import Path

from table_mae import ROOT, METHODS, PROPERTIES, calculate_maes, read_csv, number

DATASETS = ('QUEST1', 'ST_gap', 'QUEST6')
VARIANTS = (('pbe','not_applicable'), ('r2scan','not_applicable'), ('pbe0','not_applicable'),
            ('skala','on'), ('skala','off'))
FUNCTIONAL_LABELS = {'pbe':'PBE', 'r2scan':'r²SCAN', 'pbe0':'PBE0', 'skala':'Skala'}


def calculate_comparison(cp2k, pyscf, references):
    calculate_maes(cp2k, references)  # Validate full CP2K TZ coverage and gap definitions.
    cp_rows = read_csv(cp2k, ('dataset','basis','method','functional','molecule'))
    cp = {(r['dataset'],r['molecule'],r['method'],r['functional'].lower()):r
          for r in cp_rows if r['basis']=='aug-cc-pVTZ'}
    expected = {(d,m,method,xc,radii) for d,m,method,_ in cp
                for xc,radii in VARIANTS}
    seen,groups = set(),{}
    for row in read_csv(pyscf, ('dataset','basis','method','functional','molecule','radii_adjust',
                               'singlet_eV','triplet_eV','gap_eV')):
        if row['basis'] != 'aug-cc-pVTZ':
            raise ValueError('PySCF comparison supports aug-cc-pVTZ only')
        d,m,method,xc,radii = (row[k] for k in ('dataset','molecule','method','functional','radii_adjust'))
        xc=xc.lower(); key=d,m,method,xc,radii
        if key not in expected or key in seen:
            raise ValueError(f'Duplicate or unexpected PySCF result: {key}')
        seen.add(key)
        quantities = ('singlet',) if d=='QUEST6' else PROPERTIES
        values={q:number(row[q+'_eV'],str(key)) for q in quantities}
        if d=='QUEST6':
            if row['triplet_eV'] or row['gap_eV']:
                raise ValueError('QUEST6 triplet/gap fields must be empty')
        elif not math.isclose(values['gap'],values['singlet']-values['triplet'],rel_tol=0,abs_tol=5e-11):
            raise ValueError(f'Inconsistent PySCF gap: {key}')
        for q,v in values.items():
            groups.setdefault((d,method,xc,radii,q),[]).append(abs(v-float(cp[d,m,method,xc][q+'_eV'])))
    if seen != expected:
        raise ValueError(f'Missing PySCF results: {sorted(expected-seen)}')
    statistics={d:[] for d in DATASETS}
    for (d,method,xc,radii,q),values in groups.items():
        statistics[d].append(dict(method=method,functional=xc,radii_adjust=radii,quantity=q,
                                 count=len(values),mad_ev=math.fsum(values)/len(values),maxad_ev=max(values)))
    return statistics


def code_comparison_latex(statistics):
    """Render the SI summary, with PySCF Skala radii on/off and counts in the caption."""

    labels = {"QUEST1": "QUEST1", "ST_gap": "ST-gap", "QUEST6": "QUEST6"}
    counts = {}
    cells = {}
    for dataset in DATASETS:
        records = statistics[dataset]
        methods = ("UKS",) if dataset == "QUEST6" else METHODS
        quantities = ("singlet",) if dataset == "QUEST6" else PROPERTIES
        lookup = {(r["method"], r["functional"], r["radii_adjust"], r["quantity"]): r for r in records}
        expected = {(method, xc, radii, quantity) for method in methods
                    for xc, radii in VARIANTS for quantity in quantities}
        group_counts = {record["count"] for record in records}
        if len(lookup) != len(records) or set(lookup) != expected or len(group_counts) != 1:
            raise ValueError(f"Incomplete or inconsistent SI summary for {dataset}")
        counts[dataset] = group_counts.pop()
        cells[dataset] = lookup

    lines = [r"\begin{table}[htbp]", r"\centering",
             r"\caption{Mean (MAD) and maximum (MaxAD) absolute differences between PySCF",
             r"and CP2K excitation energies and singlet--triplet gaps at aug-cc-pVTZ (eV).",
             f"Statistics use {counts['QUEST1']} values per property and row for QUEST1, "
             f"{counts['ST_gap']} for ST-gap,",
             f"and {counts['QUEST6']} UKS singlet excitation energies for QUEST6.",
             r"Radii on/off refers to the PySCF Skala atomic-radius adjustment; both rows use",
             r"the same CP2K baseline. Dashes denote properties outside the reported scope.}",
             r"\label{tab:si_codes}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
             r"\begin{tabular}{@{}llrrrrrr@{}}", r"\hline",
             r"Treatment & XC functional & \multicolumn{2}{c}{$\Delta E_{\mathrm{S}}$} & "
             r"\multicolumn{2}{c}{$\Delta E_{\mathrm{T}}$} & \multicolumn{2}{c}{$\Delta E_{\mathrm{ST}}$} \\",
             r"\cline{3-4}\cline{5-6}\cline{7-8}",
             r" & & MAD & MaxAD & MAD & MaxAD & MAD & MaxAD \\"]
    for dataset in DATASETS:
        lines.extend([r"\hline", r"\multicolumn{8}{l}{\textit{" + labels[dataset] + r"}} \\"])
        methods = ("UKS",) if dataset == "QUEST6" else METHODS
        for index, method in enumerate(methods):
            if index:
                lines.append(r"\noalign{\smallskip}")
            for xc, radii in VARIANTS:
                label = FUNCTIONAL_LABELS[xc].replace("r²SCAN", r"r$^2$SCAN")
                if xc == "skala":
                    label += f" (radii {radii})"
                row = [method, label]
                for quantity in PROPERTIES:
                    if dataset == "QUEST6" and quantity != "singlet":
                        row.extend(["--", "--"])
                    else:
                        record = cells[dataset][method, xc, radii, quantity]
                        row.extend(f"{record[key]:.3f}" for key in ("mad_ev", "maxad_ev"))
                lines.append(" & ".join(row) + r" \\")
    lines.extend([r"\hline\hline", r"\end{tabular}", r"\end{table}", ""])
    return "\n".join(lines)



def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cp2k',type=Path,default=ROOT/'data/csv/cp2k_excitations.csv')
    parser.add_argument('--pyscf',type=Path,default=ROOT/'data/csv/pyscf_excitations.csv')
    parser.add_argument('--references',type=Path,default=ROOT/'data/references/excitations.csv')
    parser.add_argument('--output',type=Path,default=ROOT/'output/tables/code_comparison.tex')
    args=parser.parse_args(argv)
    try:
        cp,py,ref,output=(p.expanduser() for p in (args.cp2k,args.pyscf,args.references,args.output))
        statistics=calculate_comparison(cp,py,ref)
        comments=['% Generated by: uv run python scripts/table_codes.py',
                  '% Absolute PySCF minus CP2K differences; both Skala variants use the same CP2K baseline.']
        for label,path in [('CP2K CSV',cp),('PySCF CSV',py),('Subset references',ref)]:
            comments.append(f'% {label} SHA256: {hashlib.sha256(path.read_bytes()).hexdigest()}')
        latex='\n'.join(comments)+'\n'+code_comparison_latex(statistics)
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(latex,encoding='utf-8')
    except (ValueError,OSError,KeyError) as exc:
        parser.exit(1,f'{exc}\n')
    print('Absolute PySCF - CP2K differences (eV), aug-cc-pVTZ')
    for d in DATASETS:
        lookup={(r['method'],r['functional'],r['radii_adjust'],r['quantity']):r for r in statistics[d]}
        rows=[['Method','Functional','Radii','S MAD','S MaxAD','T MAD','T MaxAD','Gap MAD','Gap MaxAD']]
        for method in (('UKS',) if d=='QUEST6' else METHODS):
            for xc,radii in VARIANTS:
                row=[method,FUNCTIONAL_LABELS[xc],radii]
                for q in PROPERTIES:
                    r=lookup.get((method,xc,radii,q))
                    row.extend(['--','--'] if r is None else [f"{r[k]:.3f}" for k in ('mad_ev','maxad_ev')])
                rows.append(row)
        widths=[max(len(r[i]) for r in rows) for i in range(9)]
        print('\n'+d)
        for row in rows:print('  '.join(v.ljust(widths[i]) for i,v in enumerate(row)))
    print(f'\nWrote {output}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
