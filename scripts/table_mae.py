"""Print the main-text MAE table and save its LaTeX version (standard library only)."""

import argparse
import csv
import hashlib
import math
from pathlib import Path
import tempfile
import os

ROOT = Path(__file__).resolve().parents[1]
BASIS = 'aug-cc-pVTZ'
DATASETS = {'QUEST1': 'QUEST1', 'ST_gap': 'ST-gap', 'QUEST6': 'QUEST6'}
COUNTS = {'QUEST1': 12, 'ST_gap': 8, 'QUEST6': 11}
FUNCTIONALS = ('PBE', 'r2SCAN', 'PBE0', 'Skala')
METHODS = ('UKS', 'ROKS')
PROPERTIES = ('singlet', 'triplet', 'gap')


def number(value, context):
    try:
        result = float(value)
    except (ValueError, TypeError):
        raise ValueError(f'{context}: missing or invalid energy') from None
    if not math.isfinite(result):
        raise ValueError(f'{context}: nonfinite energy')
    return result


def read_csv(path, required):
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: missing required columns {required}')
        return list(reader)


def calculate_maes(calculations, references):
    refs = {}
    for row in read_csv(references, ('dataset', 'molecule', 'singlet_eV', 'triplet_eV')):
        key = row['dataset'], row['molecule']
        if key in refs or key[0] not in DATASETS:
            raise ValueError(f'Duplicate or unknown reference: {key}')
        energies = {'singlet': number(row['singlet_eV'], str(key))}
        if key[0] != 'QUEST6':
            energies['triplet'] = number(row['triplet_eV'], str(key))
            energies['gap'] = energies['singlet'] - energies['triplet']
        elif row['triplet_eV']:
            raise ValueError('QUEST6 reference triplet must be empty')
        refs[key] = energies
    for dataset, count in COUNTS.items():
        if sum(k[0] == dataset for k in refs) != count:
            raise ValueError(f'{dataset}: expected {count} reference systems')

    # A complete matched subset is required for every functional and treatment.
    expected = {(dataset, molecule, method, xc)
                for dataset, molecule in refs
                for method in (('UKS',) if dataset == 'QUEST6' else METHODS)
                for xc in FUNCTIONALS}
    seen, errors = set(), {}
    for row in read_csv(calculations, ('dataset', 'molecule', 'basis', 'method', 'functional',
                                       'singlet_eV', 'triplet_eV', 'gap_eV')):
        if row['basis'] != BASIS:
            continue
        dataset, molecule = row['dataset'], row['molecule']
        method, xc = row['method'], row['functional']
        key = dataset, molecule, method, xc
        if key in seen or key not in expected:
            raise ValueError(f'Duplicate or unexpected TZ calculation: {key}')
        seen.add(key)
        targets = refs[dataset, molecule]
        values = {p: number(row[p + '_eV'], f'{key}/{p}') for p in targets}
        if dataset == 'QUEST6':
            if row['triplet_eV'] or row['gap_eV']:
                raise ValueError(f'{key}: QUEST6 triplet/gap fields must be empty')
        elif not math.isclose(values['singlet'] - values['triplet'], values['gap'],
                              abs_tol=5e-11, rel_tol=0):
            raise ValueError(f'{key}: gap differs from singlet minus triplet')
        for prop, reference in targets.items():
            errors.setdefault((dataset, method, prop, xc), []).append(abs(values[prop] - reference))
    if seen != expected:
        raise ValueError(f'Missing TZ calculations: {sorted(expected - seen)}')
    return {key: math.fsum(values) / len(values) for key, values in errors.items()}


def table_rows(maes, latex=False):
    rows = []
    for dataset, label in DATASETS.items():
        for index, xc in enumerate(FUNCTIONALS):
            name = r'r$^2$SCAN' if latex and xc == 'r2SCAN' else xc
            cells = [label if index == 0 else '', name]
            for method in METHODS:
                for prop in PROPERTIES:
                    key = dataset, method, prop, xc
                    if key not in maes:
                        cells.append('--')
                        continue
                    value = maes[key]
                    cell = f'{value:.3f}'
                    minimum = min(maes[dataset, method, prop, f] for f in FUNCTIONALS)
                    if value == minimum:
                        cell = r'\textbf{' + cell + '}' if latex else cell + '*'
                    cells.append(cell)
            rows.append(cells)
    return rows


def latex_table(maes):
    lines = [r'''\begin{table*}[t!]
\caption{\label{tab:mae_summary}
Mean absolute errors (MAEs, in eV) of singlet ($S$) and triplet ($T$)
excitation energies and singlet--triplet gaps ($\Delta E_{\mathrm{ST}}$)
obtained with $\Delta$SCF-UKS and $\Delta$SCF-ROKS for the QUEST1, ST-gap, and QUEST6 subsets.
Open-shell singlet energies within UKS are obtained using spin purification.
The lowest MAE for each subset, property, and treatment is highlighted in bold.
For QUEST6, only UKS singlet excitation-energy MAEs are reported.
}
\begin{ruledtabular}
\begin{tabular}{llcccccc}
 & & \multicolumn{3}{c}{$\Delta$SCF-UKS} & \multicolumn{3}{c}{$\Delta$SCF-ROKS} \\
\cline{3-5}\cline{6-8}
Dataset & XC functional & $S$ & $T$ & $\Delta E_{\mathrm{ST}}$ & $S$ & $T$ & $\Delta E_{\mathrm{ST}}$ \\
\hline''']
    for i, row in enumerate(table_rows(maes, latex=True)):
        if i and i % len(FUNCTIONALS) == 0:
            lines.append(r'\hline')
        lines.append(' & '.join(row) + r' \\')
    lines += [r'\end{tabular}', r'\end{ruledtabular}', r'\end{table*}', '']
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculations', type=Path, default=ROOT / 'data/csv/cp2k_excitations.csv')
    parser.add_argument('--references', type=Path, default=ROOT / 'data/references/excitations.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/tables/mae_summary.tex')
    args = parser.parse_args(argv)
    try:
        calc, ref, output = (p.expanduser() for p in (args.calculations, args.references, args.output))
        maes = calculate_maes(calc, ref)
        comments = ['% Generated by: uv run python scripts/table_mae.py',
                    f'% CP2K / {BASIS}; full subsets; minima selected before rounding.']
        for label, path in [('Calculation CSV', calc), ('Reference CSV', ref)]:
            comments.append(f'% {label} SHA256: {hashlib.sha256(path.read_bytes()).hexdigest()}')
        content = '\n'.join(comments) + '\n' + latex_table(maes)
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
    except (OSError, ValueError) as exc:
        parser.exit(1, f'{exc}\n')
    rows = [['Dataset', 'Functional', 'UKS S', 'UKS T', 'UKS gap', 'ROKS S', 'ROKS T', 'ROKS gap'],
            *table_rows(maes)]
    widths = [max(len(row[i]) for row in rows) for i in range(8)]
    print(f'MAEs (eV), CP2K / {BASIS}\n')
    for row in rows:
        print('  '.join(cell.ljust(widths[i]) for i, cell in enumerate(row)))
    print('\n* Lowest unrounded MAE within each dataset, property and treatment.')
    print(f'\nWrote {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
