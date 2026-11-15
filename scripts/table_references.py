"""Generate the three SM reference tables from reference CSV and state metadata.

Uses only Python's standard library; no development-project dependencies.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {'QUEST1': 12, 'ST_gap': 8, 'QUEST6': 11}

def selection_tables_latex(records):
    """Render reference labels separately from numerical calculation results."""
    symmetries = {"A''": r"A^{\prime\prime}", "A1": "A_1", "A2": "A_2",
                  "B1": "B_1", "B2": "B_2", "B1u": "B_{1u}"}
    characters = {"n->pi*": r"$n\to\pi^*$", "pi->pi*": r"$\pi\to\pi^*$",
                  "n->3s": r"$n\to3s$"}

    def state(symmetry, multiplicity):
        return rf"${{}}^{{{multiplicity}}}{symmetries[symmetry]}$"

    def table(caption, label, columns, header, rows):
        lines = [r"\begin{table}[htbp]", r"\centering", r"\caption{" + caption + "}",
                 r"\label{" + label + "}", r"\small", r"\setlength{\tabcolsep}{6pt}",
                 r"\begin{tabular}{@{}" + columns + "@{}}", r"\hline", header + r" \\",
                 r"\hline"]
        lines.extend(" & ".join(row) + r" \\" for row in rows)
        lines.extend([r"\hline\hline", r"\end{tabular}", r"\end{table}", ""])
        return "\n".join(lines)

    quest1 = []
    for entry in records["QUEST1"]:
        kind = {"valence": "V", "Rydberg": "R"}[entry["kind"]]
        name = entry["system"]
        if entry["charge"]:
            name += rf" ($q={entry['charge']:+d}$)"
        quest1.append([name, state(entry["symmetry"], 1), state(entry["symmetry"], 3),
                       kind + "; " + characters[entry["singlet_character"]],
                       f"{entry['tbe_ev']['singlet']:.2f}", f"{entry['tbe_ev']['triplet']:.2f}"])
    tables = {"selection_quest1.tex": table(
        r"Selected QUEST1 targets and corrected TBEs (eV) from Table~6 of Ref.~\onlinecite{Loos2018}. "
        r"Character refers to the singlet; V and R denote valence and Rydberg excitations.",
        "tab:si_states_quest1", "lcccrr",
        r"System & Singlet & Triplet & Singlet character & $\Delta E_{\mathrm{S}}^{\mathrm{TBE}}$ & $\Delta E_{\mathrm{T}}^{\mathrm{TBE}}$",
        quest1)}

    st_gap = []
    for entry in records["ST_gap"]:
        number = entry["compound"]
        name = (f"{entry['system']} ({number})" if number <= 2 else f"Molecule~{number}")
        st_gap.append([name, f"{entry['tbe_ev']['singlet']:.3f}",
                       f"{entry['tbe_ev']['triplet']:.3f}", f"{entry['tbe_ev']['gap']:.3f}"])
    tables["selection_st_gap.tex"] = table(
        r"Selected ST-gap systems and vertical TBEs (eV) for $S_1$ and $T_1$ "
        r"from Table~S1 in the SM of Ref.~\onlinecite{Loos2023}. "
        r"Compound numbering follows that reference; $\Delta E_{\mathrm{ST}}=\Delta E_{\mathrm{S}}-\Delta E_{\mathrm{T}}$.",
        "tab:si_states_st_gap", "lrrr",
        r"System & $\Delta E_{\mathrm{S}}^{\mathrm{TBE}}$ & $\Delta E_{\mathrm{T}}^{\mathrm{TBE}}$ & $\Delta E_{\mathrm{ST}}^{\mathrm{TBE}}$",
        st_gap)

    quest6 = []
    for entry in records["QUEST6"]:
        target = rf"$S_{{{entry['singlet_root']}}}$ (" + state(entry["symmetry"], 1) + ")"
        name = entry["system"]
        if name == "Dimethylaminobenzonitrile":
            name += " (DMABN)"
        elif name == "N-Phenylpyrrole":
            name = r"$N$-Phenylpyrrole (PP)"
        elif name == "Nitropyridine N-Oxide":
            name = r"Nitropyridine $N$-oxide"
        quest6.append([name, target, characters[entry["singlet_character"]],
                       f"{entry['tbe_ev']['singlet']:.2f}"])
    tables["selection_quest6.tex"] = table(
        r"Selected QUEST6 singlet CT targets and aug-cc-pVQZ TBEs (eV) from "
        r"Table~1 of Ref.~\onlinecite{Loos2021}. The $S_n$ labels follow the "
        r"LR-CCSD/cc-pVTZ ordering in its SM (Table~S1).",
        "tab:si_states_quest6", "lccr",
        r"System & Reference singlet & Character & $\Delta E_{\mathrm{S}}^{\mathrm{TBE}}$",
        quest6)
    return tables



def load_references(energies_path, states_path):
    metadata = json.loads(states_path.read_text())['datasets']
    if set(metadata) != set(COUNTS):
        raise ValueError('State metadata must cover QUEST1, ST_gap and QUEST6')
    energies = {}
    with energies_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        required = {'dataset', 'molecule', 'singlet_eV', 'triplet_eV', 'source'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Missing reference CSV columns')
        for row in reader:
            key = row['dataset'], row['molecule']
            if key in energies or key[0] not in COUNTS:
                raise ValueError(f'Duplicate or unknown reference: {key}')
            values = {}
            for prop in ('singlet', 'triplet'):
                value = row[prop + '_eV']
                if prop == 'triplet' and key[0] == 'QUEST6':
                    if value:
                        raise ValueError('QUEST6 must not contain reference triplets')
                    values[prop] = None
                else:
                    values[prop] = float(value)
                    if not math.isfinite(values[prop]):
                        raise ValueError(f'Nonfinite reference: {key}/{prop}')
            values['gap'] = None if values['triplet'] is None else values['singlet'] - values['triplet']
            if row['source'] != 'https://doi.org/' + metadata[key[0]]['doi']:
                raise ValueError(f'Reference source mismatch: {key}')
            energies[key] = values
    records, seen = {}, set()
    for dataset, count in COUNTS.items():
        entries = metadata[dataset]['entries']
        if len(entries) != count:
            raise ValueError(f'{dataset}: expected {count} states')
        records[dataset] = []
        for entry in entries:
            key = dataset, entry['molecule']
            if key in seen or key not in energies:
                raise ValueError(f'Duplicate state or missing energy: {key}')
            if entry['system'].replace(' ', '_') != entry['molecule']:
                raise ValueError(f'Inconsistent molecule identifier: {key}')
            seen.add(key)
            records[dataset].append(dict(entry, tbe_ev=energies[key]))
    if seen != set(energies):
        raise ValueError('Reference energies without matching state metadata')
    # The formatter uses the manuscript citation keys; reject mismatched metadata.
    expected_keys = {'QUEST1':'Loos2018','ST_gap':'Loos2023','QUEST6':'Loos2021'}
    if any(metadata[d]['bibkey'] != key for d, key in expected_keys.items()):
        raise ValueError('State metadata citation keys differ from table captions')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references', type=Path, default=ROOT / 'data/references/excitations.csv')
    parser.add_argument('--states', type=Path, default=ROOT / 'data/references/states.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'output/tables')
    parser.add_argument('--print', action='store_true', dest='print_tables', help='Also print the generated LaTeX')
    args = parser.parse_args(argv)
    try:
        energies, states = args.references.expanduser(), args.states.expanduser()
        records = load_references(energies, states)
        tables = selection_tables_latex(records)
        header = '% Generated by: uv run python scripts/table_references.py\n'
        for label, path in [('Reference CSV', energies), ('State metadata', states)]:
            header += f'% {label} SHA256: {hashlib.sha256(path.read_bytes()).hexdigest()}\n'
        output = args.output_dir.expanduser()
        output.mkdir(parents=True, exist_ok=True)
        for name, table in tables.items():
            (output / name).write_text(header + table, encoding='utf-8')
            print(f'Wrote {output / name}')
            if args.print_tables:
                print(table)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f'{exc}\n')
    print('Validated 31 reference systems and 51 excitation TBEs.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
