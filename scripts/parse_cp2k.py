"""Export CP2K excitation energies using only the Python standard library.

Input layout: data/cp2k/<basis>/<dataset>/<uks|roks>/<functional>/<molecule>/
Each molecular directory must contain ground_state.out, singlet.out and triplet.out.
UKS: E(S) = 2 E(mixed) - E(triplet). ROKS: use the singlet energy directly.
The triplet energy always comes from the separate triplet.out calculation.
Only benchmarked singlets are exported for QUEST6, but its triplet is still
required for UKS spin purification. No references, exclusions or MAEs are applied.
Numerical SCF completion does not establish spectroscopic state assignment.
"""

import argparse
from collections import Counter
import csv
import math
import os
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[1]
# Conversion used for the published benchmark exports (Ha -> eV).
HARTREE_TO_EV = 27.211396132
STAGES = ("ground_state", "singlet", "triplet")
FIELDS = ("dataset", "molecule", "basis", "method", "functional",
          "singlet_eV", "triplet_eV", "gap_eV")
FUNCTIONALS = {"pbe": "PBE", "pbe0": "PBE0", "r2scan": "r2SCAN", "skala": "Skala"}
ENERGY = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*QS\s*\) energy \[hartree\]\s+(\S+)")
CONVERGED = re.compile(r"(?:SCF run|outer SCF loop)\s+converged\b")
FAILED = re.compile(r"(?:SCF run|outer SCF loop)\s+NOT converged\b", re.I)
SPIN = re.compile(r"^\s*DFT\|.*\b(ROKS|UKS|RKS)\s*$", re.M)


def parse_output(path, expected_spin):
    """Read the latest invocation only; require converged SCF and normal exit."""
    text = path.read_text(encoding="utf-8", errors="replace").rsplit("PROGRAM STARTED AT", 1)[-1]
    if FAILED.search(text) or "[ABORT]" in text or "MPI_ABORT" in text.upper():
        raise ValueError("SCF failure or aborted calculation")
    spin = set(SPIN.findall(text))
    if spin != {expected_spin}:
        raise ValueError(f"expected {expected_spin} header, found {sorted(spin)}")
    converged = list(CONVERGED.finditer(text))
    if not converged or converged[-1].start() < text.rfind("SCF WAVEFUNCTION OPTIMIZATION"):
        raise ValueError("no final SCF convergence marker")
    energies = list(ENERGY.finditer(text))
    if not energies or energies[-1].start() < converged[-1].end():
        raise ValueError("no total energy after final SCF convergence")
    energy = float(energies[-1].group(1).replace("D", "E").replace("d", "e"))
    if not math.isfinite(energy):
        raise ValueError("nonfinite total energy")
    if text.rfind("PROGRAM ENDED AT") < energies[-1].end():
        raise ValueError("no normal termination after final energy")
    return energy


def collect(data_dir):
    """Parse all molecular directories; fail rather than silently omit bad data."""
    if not data_dir.is_dir():
        raise ValueError(f"input directory not found: {data_dir}")
    directories = sorted(p for p in data_dir.glob("*/*/*/*/*") if p.is_dir())
    if not directories:
        raise ValueError(f"no molecular directories found under {data_dir}")
    rows, errors = [], []
    for directory in directories:
        basis, dataset, method, functional, molecule = directory.relative_to(data_dir).parts
        if (dataset not in {"QUEST1", "ST_gap", "QUEST6"} or method not in {"uks", "roks"}
                or functional not in FUNCTIONALS):
            errors.append(f"{directory}: unsupported dataset, method or functional")
            continue
        if dataset == "QUEST6" and method != "uks":
            errors.append(f"{directory}: QUEST6 exports are defined for UKS only")
            continue
        energies = {}
        for stage in STAGES:
            path = directory / f"{stage}.out"
            try:
                energies[stage] = parse_output(path, "RKS" if stage == "ground_state" else method.upper())
            except (OSError, ValueError) as exc:
                errors.append(f"{path}: {exc}")
        if len(energies) != len(STAGES):
            continue
        ground, singlet, triplet = (energies[s] for s in STAGES)
        if method == "uks":
            singlet = 2 * singlet - triplet
        singlet_ev = (singlet - ground) * HARTREE_TO_EV
        triplet_ev = (triplet - ground) * HARTREE_TO_EV
        rows.append(dict(zip(FIELDS, (
            dataset, molecule, basis, method.upper(), FUNCTIONALS[functional],
            singlet_ev, triplet_ev if dataset != "QUEST6" else "",
            singlet_ev - triplet_ev if dataset != "QUEST6" else "",
        ))))
    if errors:
        raise ValueError("No CSV written. Invalid or incomplete calculations:\n" + "\n".join(errors))
    return sorted(rows, key=lambda r: tuple(r[k] for k in FIELDS[:5]))


def write_csv(rows, output):
    """Replace the CSV only after all data have been validated."""
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8",
                                         dir=output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: format(v, ".12f") if isinstance(v, float) else v
                                 for k, v in row.items()})
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/cp2k")
    parser.add_argument("--output", type=Path, default=ROOT / "data/csv/cp2k_excitations.csv")
    args = parser.parse_args(argv)
    try:
        rows = collect(args.data_dir.expanduser().resolve())
        output = args.output.expanduser().resolve()
        write_csv(rows, output)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{exc}\n")
    print(f"Wrote {len(rows)} rows to {output}")
    for key, count in sorted(Counter((r['dataset'], r['basis'], r['method']) for r in rows).items()):
        print(f"  {' / '.join(key)}: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
