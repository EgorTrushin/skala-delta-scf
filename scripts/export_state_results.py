"""Join benchmark energies to TBEs and signed errors, without rounding."""

import argparse
import csv
from pathlib import Path

from table_mae import ROOT, PROPERTIES, read_csv, number
from table_codes import calculate_comparison as validate_codes
from table_basis import calculate_comparison as validate_bases


def collect_results(cp2k, pyscf, references):
    validate_codes(cp2k, pyscf, references)
    validate_bases(cp2k, references)
    refs = {(r['dataset'], r['molecule']): r for r in
            read_csv(references, ('dataset', 'molecule', 'source'))}
    result = []
    for code, path in [('CP2K', cp2k), ('PySCF', pyscf)]:
        for row in read_csv(path, ('dataset', 'molecule', 'basis')):
            if row['basis'] not in ('aug-cc-pVTZ', 'aug-cc-pVQZ') or (
                    row['dataset'] == 'QUEST6' and row['basis'] != 'aug-cc-pVTZ'):
                raise ValueError(f'Unsupported basis/dataset: {row}')
            ref = refs[row['dataset'], row['molecule']]
            targets = {'singlet': number(ref['singlet_eV'], str(ref))}
            if row['dataset'] != 'QUEST6':
                targets['triplet'] = number(ref['triplet_eV'], str(ref))
                targets['gap'] = targets['singlet'] - targets['triplet']
            for prop in PROPERTIES:
                if prop not in targets:
                    continue
                value = number(row[prop + '_eV'], str(row))
                result.append(dict(code=code, **{k: row[k] for k in
                    ('dataset', 'molecule', 'basis', 'method', 'functional')},
                    radii_adjust=row.get('radii_adjust', 'not_applicable'),
                    property=prop, calculated_eV=value, tbe_eV=targets[prop],
                    error_eV=value-targets[prop], reference_source=ref['source']))
    return result


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cp2k', type=Path, default=ROOT/'data/csv/cp2k_excitations.csv')
    parser.add_argument('--pyscf', type=Path, default=ROOT/'data/csv/pyscf_excitations.csv')
    parser.add_argument('--references', type=Path, default=ROOT/'data/references/excitations.csv')
    parser.add_argument('--output', type=Path, default=ROOT/'data/csv/state_results.csv')
    args = parser.parse_args()
    try:
        rows = collect_results(args.cp2k, args.pyscf, args.references)
        write_csv(args.output, rows)
    except (OSError, ValueError) as exc:
        parser.exit(1, f'{exc}\n')
    print(f'Wrote {len(rows)} energies/gaps with TBEs and signed errors to {args.output}')


if __name__ == '__main__':
    main()
