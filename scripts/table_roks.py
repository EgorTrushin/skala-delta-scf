"""Print and export the two ROKS/ST-gap supplementary tables from the validated CSV."""
import argparse
import csv
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ['Heptazine', 'Cyclazine', *[f'Molecule{i}' for i in range(3, 9)]]
FUNCTIONALS = ['pbe', 'r2scan', 'pbe0', 'skala']
LABELS = {'pbe': 'PBE', 'r2scan': r'r$^2$SCAN', 'pbe0': 'PBE0', 'skala': 'Skala'}
PAIRS = [('pbe', 'pbe'), ('skala', 'pbe'), ('pbe', 'skala'), ('skala', 'skala')]
PROPERTIES = ['common_gap_eV', 'triplet_correction_eV', 'gap_eV']


def load_data(path):
    expected = {(m, e, d) for m in MOLECULES for e, d in PAIRS + [('r2scan', 'r2scan'), ('pbe0', 'pbe0')]}
    data = {}
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            key = row['molecule'], row['evaluating_functional'], row['orbital_functional']
            if key not in expected or key in data or row['basis'] != 'aug-cc-pVTZ':
                raise ValueError(f'Duplicate or unexpected row: {key}')
            values = {p: float(row[p]) for p in PROPERTIES}
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError(f'Nonfinite energy: {key}')
            if abs(values['common_gap_eV'] + values['triplet_correction_eV'] - values['gap_eV']) > 1e-10:
                raise ValueError(f'Inconsistent gap: {key}')
            data[key] = values
    if set(data) != expected:
        raise ValueError(f'Incomplete input: missing {sorted(expected - set(data))}')
    return data


def table_rows(data, columns):
    rows = [(m, [data[m, e, d][p] for e, d, p in columns]) for m in MOLECULES]
    means = [math.fsum(row[1][i] for row in rows) / len(rows) for i in range(len(columns))]
    return rows + [('Mean', means)]


def value(x, latex=False):
    text = f'{x:+.3f}'
    if text == '-0.000':
        text = '+0.000'
    return text.replace('-', r'\textminus{}') if latex else text


def body(rows):
    lines = []
    for molecule, values in rows:
        if molecule == 'Mean':
            lines.append(r'\hline')
        name = molecule.replace('Molecule', 'Molecule~')
        lines.append(' & '.join([name, *[value(v, True) for v in values]]) + r' \\')
    return lines


def inverted(data, evaluator, donor, prop):
    return sum(data[m, evaluator, donor][prop] < 0 for m in MOLECULES)


def render_frozen(data):
    columns = [(e, d, 'gap_eV') for e, d in PAIRS]
    rows = table_rows(data, columns)
    counts = ', '.join(str(inverted(data, e, d, 'gap_eV')) for e, d in PAIRS)
    caption = (r'Frozen-orbital ROKS singlet--triplet gaps for ST-gap (CP2K, aug-cc-pVTZ; eV). '
               r'$F@G$ denotes evaluation with $F$ on the respective singlet and separately '
               r'optimized restricted-triplet orbitals obtained with $G$. '
               r'The four columns contain ' + counts + r' negative gaps, respectively, among the eight molecules. '
               r'Both functional/orbital comparison paths are retained; their contributions are path-dependent.')
    lines = [r'\begin{table}[htbp]', r'\centering', r'\caption{' + caption + '}',
             r'\label{tab:si_roks_frozen_orbitals}', r'\footnotesize',
             r'\setlength{\tabcolsep}{5pt}', r'\begin{tabular}{@{}lrrrr@{}}', r'\hline\hline',
             r'Molecule & PBE@PBE & Skala@PBE & PBE@Skala & Skala@Skala \\', r'\hline']
    return '\n'.join(lines + body(rows) + [r'\hline\hline', r'\end{tabular}', r'\end{table}', '']), rows


def render_triplet(data):
    columns = [(f, f, p) for f in FUNCTIONALS for p in ('common_gap_eV', 'gap_eV')]
    rows = table_rows(data, columns)
    common = inverted(data, 'skala', 'skala', 'common_gap_eV')
    optimized = inverted(data, 'skala', 'skala', 'gap_eV')
    conventional = any(inverted(data, f, f, p) for f in FUNCTIONALS[:-1] for p in ('common_gap_eV', 'gap_eV'))
    inversion_text = (f'Skala gives {common} common-orbital and {optimized} optimized-gap inversions among the eight molecules.'
                      + (' No inversions occur for the other functionals.' if not conventional else ''))
    caption = (r'Effect of separate triplet optimization on ROKS singlet--triplet gaps for ST-gap '
               r'(CP2K, aug-cc-pVTZ; eV). Common: $E_{\mathrm S}(C_{\mathrm S})-E_{\mathrm T}(C_{\mathrm S})$; '
               r'optimized: $E_{\mathrm S}(C_{\mathrm S})-E_{\mathrm T}(C_{\mathrm T})$, '
               r'where $C_{\mathrm S}$ and $C_{\mathrm T}$ are singlet- and restricted-triplet-optimized orbitals. '
               r'The singlet and nuclear geometry remain fixed. ' + inversion_text)
    lines = [r'\begin{table}[htbp]', r'\centering', r'\caption{' + caption + '}',
             r'\label{tab:si_triplet_relaxation}', r'\footnotesize',
             r'\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{@{}lrrrrrrrr@{}}', r'\hline\hline',
             'Molecule & ' + ' & '.join(r'\multicolumn{2}{c}{' + LABELS[f] + '}' for f in FUNCTIONALS) + r' \\',
             r'\cline{2-3}\cline{4-5}\cline{6-7}\cline{8-9}',
             ' & ' + ' & '.join(['Common', 'Optimized'] * 4) + r' \\', r'\hline']
    return '\n'.join(lines + body(rows) + [r'\hline\hline', r'\end{tabular}', r'\end{table}', '']), rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data/csv/roks_st_gap.csv')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'output/tables')
    args = parser.parse_args()
    data = load_data(args.data)
    header = ('% Generated by: uv run python scripts/table_roks.py\n'
              '% Source CSV SHA256: ' + hashlib.sha256(args.data.read_bytes()).hexdigest() + '\n'
              '% Means and inversion counts use unrounded energies.\n')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, render in [('roks_frozen_orbitals.tex', render_frozen),
                             ('roks_triplet_relaxation.tex', render_triplet)]:
        latex, rows = render(data)
        path = args.output_dir / filename
        path.write_text(header + latex)
        print('\n' + filename + ' (eV)')
        for name, values in rows:
            print(f'{name:12s} ' + '  '.join(value(v) for v in values))
        print(f'Wrote {path}')

if __name__ == '__main__':
    main()
