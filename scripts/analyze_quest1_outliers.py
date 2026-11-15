"""Reproduce the descriptive QUEST1 sensitivity analysis (CP2K/aug-cc-pVTZ).

The three named exclusions are applied identically to all comparisons. They
are not a statistical rejection rule; the full benchmark remains primary.
"""

import argparse
import math
from pathlib import Path

from table_mae import ROOT, BASIS, METHODS, FUNCTIONALS, PROPERTIES, calculate_maes, read_csv
from export_state_results import write_csv

EXCLUDED = ('Ethylene', 'Cyclopropene', 'Streptocyanine-C1')


def analyze(calculations, references):
    maes = calculate_maes(calculations, references)
    refs = {r['molecule']: r for r in read_csv(references, ('dataset', 'molecule'))
            if r['dataset'] == 'QUEST1'}
    if not set(EXCLUDED) <= refs.keys():
        raise ValueError('Missing named exclusion in reference data')
    groups, details = {}, []
    for row in read_csv(calculations, ('dataset', 'basis')):
        if (row['dataset'], row['basis']) != ('QUEST1', BASIS):
            continue
        ref = refs[row['molecule']]
        targets = [float(ref['singlet_eV']), float(ref['triplet_eV'])]
        targets.append(targets[0] - targets[1])
        for prop, tbe in zip(PROPERTIES, targets):
            error = float(row[prop+'_eV']) - tbe
            excluded = row['molecule'] in EXCLUDED
            groups.setdefault((row['method'], row['functional'], prop), []).append((excluded, error))
            details.append(dict(molecule=row['molecule'], method=row['method'],
                                functional=row['functional'], property=prop,
                                excluded=excluded, calculated_eV=float(row[prop+'_eV']),
                                tbe_eV=tbe, error_eV=error))
    summary = []
    for method in METHODS:
        for xc in FUNCTIONALS:
            for prop in PROPERTIES:
                values = groups[method, xc, prop]
                retained = [abs(e) for excluded, e in values if not excluded]
                removed = math.fsum(abs(e) for excluded, e in values if excluded)
                total = math.fsum(abs(e) for _, e in values)
                if len(values) != 12 or len(retained) != 9:
                    raise ValueError('Unmatched QUEST1 subset')
                summary.append(dict(method=method, functional=xc, property=prop,
                    full_count=12, retained_count=9,
                    full_mae_eV=maes['QUEST1', method, prop, xc],
                    retained_mae_eV=math.fsum(retained)/9,
                    excluded_absolute_error_eV=removed,
                    excluded_fraction_absolute_error=removed/total if total else 0.0))
    return summary, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculations', type=Path, default=ROOT/'data/csv/cp2k_excitations.csv')
    parser.add_argument('--references', type=Path, default=ROOT/'data/references/excitations.csv')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'output/tables')
    args = parser.parse_args()
    try:
        summary, details = analyze(args.calculations, args.references)
        write_csv(args.output_dir/'quest1_outlier_summary.csv', summary)
        write_csv(args.output_dir/'quest1_outlier_errors.csv', details)
    except (OSError, ValueError) as exc:
        parser.exit(1, f'{exc}\n')
    print('Common exclusions: ' + ', '.join(EXCLUDED))
    print('MAEs in eV; full 12-molecule benchmark remains primary.\n')
    print('Method Functional Property   Full (12) Retained (9)')
    for row in summary:
        print(f"{row['method']:6} {row['functional']:10} {row['property']:8} "
              f"{row['full_mae_eV']:9.3f} {row['retained_mae_eV']:12.3f}")
        if (row['method'], row['functional'], row['property']) == ('ROKS', 'Skala', 'gap'):
            print(f"  Excluded share of absolute gap error: {100*row['excluded_fraction_absolute_error']:.2f}%")
    print(f'\nWrote summary and molecular errors to {args.output_dir}')


if __name__ == '__main__':
    main()
