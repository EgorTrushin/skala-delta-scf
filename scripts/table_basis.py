"""Print and export the CP2K TZ/QZ basis comparison (Python standard library only)."""

import argparse
import hashlib
import math
import os
from pathlib import Path
import tempfile

from table_mae import ROOT, FUNCTIONALS, METHODS, PROPERTIES, read_csv, number

BASES = ('aug-cc-pVTZ', 'aug-cc-pVQZ')
DATASETS = {'QUEST1': 'QUEST1', 'ST_gap': 'ST-gap'}
COUNTS = {'QUEST1': 12, 'ST_gap': 8}


def calculate_comparison(calculations, references):
    # References identify the complete subsets; no reference energies enter MADs.
    systems = {dataset: set() for dataset in DATASETS}
    for row in read_csv(references, ('dataset', 'molecule')):
        dataset, molecule = row['dataset'], row['molecule']
        if dataset not in systems:
            continue
        if molecule in systems[dataset]:
            raise ValueError(f'Duplicate reference system: {dataset}/{molecule}')
        systems[dataset].add(molecule)
    for dataset, count in COUNTS.items():
        if len(systems[dataset]) != count:
            raise ValueError(f'{dataset}: expected {count} reference systems')
    cells = {}
    for row in read_csv(calculations, ('dataset', 'molecule', 'basis', 'method', 'functional',
                                       'singlet_eV', 'triplet_eV', 'gap_eV')):
        dataset = row['dataset']
        if dataset not in DATASETS:
            continue
        basis, method, xc, molecule = (row[k] for k in ('basis', 'method', 'functional', 'molecule'))
        if basis not in BASES:
            continue
        key = dataset, basis, method, xc, molecule
        if (method not in METHODS or xc not in FUNCTIONALS or molecule not in systems[dataset]
                or key in cells):
            raise ValueError(f'Duplicate or unexpected calculation: {key}')
        values = {p: number(row[p + '_eV'], f'{key}/{p}') for p in PROPERTIES}
        if not math.isclose(values['singlet'] - values['triplet'], values['gap'],
                            rel_tol=0, abs_tol=5e-11):
            raise ValueError(f'Inconsistent gap: {key}')
        cells[key] = values
    statistics, missing = {}, []
    for dataset in DATASETS:
        for method in METHODS:
            for xc in FUNCTIONALS:
                groups = []
                for basis in BASES:
                    group = {key[-1]: value for key, value in cells.items()
                             if key[:4] == (dataset, basis, method, xc)}
                    if not group and basis == BASES[1] and (dataset, xc) == ('ST_gap', 'PBE0'):
                        missing.append((dataset, method, xc))
                    elif set(group) != systems[dataset]:
                        raise ValueError(f'Incomplete {dataset}/{basis}/{method}/{xc}: '
                                         f'missing {sorted(systems[dataset] - set(group))}')
                    groups.append(group)
                if not groups[1]:
                    continue
                for prop in PROPERTIES:
                    differences = [abs(groups[1][m][prop] - groups[0][m][prop])
                                   for m in sorted(systems[dataset])]
                    statistics[dataset, method, xc, prop] = (
                        math.fsum(differences) / len(differences), max(differences))
    return statistics, missing


def rows_for_dataset(statistics, dataset, latex=False):
    rows = []
    for method in METHODS:
        for xc in FUNCTIONALS:
            row = [method, r'r$^2$SCAN' if latex and xc == 'r2SCAN' else xc]
            for prop in PROPERTIES:
                values = statistics.get((dataset, method, xc, prop))
                row.extend(['--', '--'] if values is None else [f'{v:.3f}' for v in values])
            rows.append(row)
    return rows


def render_table(statistics, missing):
    caption = ('Mean absolute differences (MADs) and maximum absolute differences (MaxADs) '
               'between CP2K aug-cc-pVTZ and aug-cc-pVQZ excitation energies and '
               'singlet--triplet gaps (eV). Statistics include 12 molecules for QUEST1 '
               'and eight for ST-gap, where available. ')
    if missing:
        methods = {m for _, m, _ in missing}
        qualifier = '' if methods == set(METHODS) else ' (' + ', '.join(sorted(methods)) + ')'
        caption += 'Dashes indicate unavailable aug-cc-pVQZ results for ST-gap PBE0' + qualifier + '. '
    caption += 'No aug-cc-pVQZ calculations were performed for QUEST6.'
    lines = [r'\begin{table}[htbp]', r'\centering', r'\caption{' + caption + '}',
             r'\label{tab:si_basis}', r'\footnotesize', r'\setlength{\tabcolsep}{4pt}',
             r'\begin{tabular}{@{}llrrrrrr@{}}', r'\hline',
             r'Treatment & XC functional & \multicolumn{2}{c}{$\Delta E_{\mathrm{S}}$} & \multicolumn{2}{c}{$\Delta E_{\mathrm{T}}$} & \multicolumn{2}{c}{$\Delta E_{\mathrm{ST}}$} \\',
             r'\cline{3-4}\cline{5-6}\cline{7-8}',
             r' & & MAD & MaxAD & MAD & MaxAD & MAD & MaxAD \\', r'\hline']
    for index, (dataset, label) in enumerate(DATASETS.items()):
        if index:
            lines.append(r'\hline')
        lines.append(r'\multicolumn{8}{l}{\textit{' + label + r'}} \\')
        for i, row in enumerate(rows_for_dataset(statistics, dataset, latex=True)):
            if i == len(FUNCTIONALS):
                lines.append(r'\noalign{\smallskip}')
            lines.append(' & '.join(row) + r' \\')
    return '\n'.join(lines + [r'\hline\hline', r'\end{tabular}', r'\end{table}', ''])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculations', type=Path, default=ROOT / 'data/csv/cp2k_excitations.csv')
    parser.add_argument('--references', type=Path, default=ROOT / 'data/references/excitations.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/tables/basis_comparison.tex')
    args = parser.parse_args(argv)
    try:
        calc, ref, output = (p.expanduser() for p in (args.calculations, args.references, args.output))
        statistics, missing = calculate_comparison(calc, ref)
        comments = ['% Generated by: uv run python scripts/table_basis.py',
                    '% Absolute paired QZ minus TZ differences; no TBE errors or rounded inputs.']
        for label, path in [('Calculation CSV', calc), ('Subset reference CSV', ref)]:
            comments.append(f'% {label} SHA256: {hashlib.sha256(path.read_bytes()).hexdigest()}')
        content = '\n'.join(comments) + '\n' + render_table(statistics, missing)
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=output.parent,
                                             delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(content)
            os.replace(temporary, output)
        finally:
            if temporary:
                temporary.unlink(missing_ok=True)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'{exc}\n')
    print('CP2K basis comparison: absolute QZ - TZ differences (eV)')
    for dataset, label in DATASETS.items():
        rows = [['Method', 'Functional', 'S MAD', 'S MaxAD', 'T MAD', 'T MaxAD', 'Gap MAD', 'Gap MaxAD'],
                *rows_for_dataset(statistics, dataset)]
        widths = [max(len(row[i]) for row in rows) for i in range(8)]
        print('\n' + label)
        for row in rows:
            print('  '.join(cell.ljust(widths[i]) for i, cell in enumerate(row)))
    if missing:
        print('\nUnavailable QZ groups: ' + ', '.join('/'.join(g) for g in missing))
    print(f'\nWrote {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
