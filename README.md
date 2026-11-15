# Skala ΔSCF benchmark data

Data and scripts for reproducing the figures and tables in
*Beyond the Ground State: Assessing the Machine-Learned Skala Functional in
ΔSCF Calculations of Molecular Excitation Energies* and its Supplementary Material.

Run all commands from the repository root:

```bash
uv sync --locked
```

The scripts use the supplied CSV files in `data/csv/` and reference data in
`data/references/`. CP2K and PySCF installations are not needed to reproduce
these figures and tables.

## Figures

```bash
uv run python scripts/plot_quest1.py
uv run python scripts/plot_st_gap.py
uv run python scripts/plot_st_gap_frozen.py
uv run python scripts/plot_quest6.py
```

PDFs are saved in `output/figures/`:

| Script | Output |
| --- | --- |
| `plot_quest1.py` | `quest1_signed_errors.pdf` |
| `plot_st_gap.py` | `st_gap_signed_errors.pdf`, `st_gap_gap_distributions.pdf` |
| `plot_st_gap_frozen.py` | `st_gap_roks_analysis.pdf` |
| `plot_quest6.py` | `quest6_scatter.pdf` |

Add `--png` for PNG previews. The frozen-orbital figure uses
`data/csv/roks_st_gap.csv` and preserves the original PBE–Skala donor snapshots.

## Tables

```bash
uv run python scripts/table_mae.py          # Main-text MAEs
uv run python scripts/table_references.py   # SM: selected states and TBEs
uv run python scripts/table_basis.py        # SM: aug-cc-pVTZ / aug-cc-pVQZ
uv run python scripts/table_codes.py        # SM: CP2K / PySCF
uv run python scripts/table_roks.py         # SM: frozen orbitals and triplet relaxation
```

LaTeX files are saved in `output/tables/`. The scripts also print numerical
tables to the terminal; use `--print` with `table_references.py` to print its
LaTeX. Each ROKS table is exported as a separate file.

For the additional comparison with published PBE excitation energies:

```bash
uv run python scripts/table_pbe_literature.py
```

This prints the molecular comparisons and saves their CSV tables in
`output/tables/`. All scripts provide `--help` for alternative input/output paths.

For molecular energies, TBEs and signed errors, and the QUEST1 sensitivity
analysis discussed in the main text:

```bash
uv run python scripts/export_state_results.py
uv run python scripts/analyze_quest1_outliers.py
```

The first writes `data/csv/state_results.csv`: one row per energy or gap, all
available basis sets, both PySCF Skala radii settings, and reference sources.
Errors are calculated minus TBE; reference gaps are singlet minus triplet TBEs.
QUEST6 contains singlets only. Frozen probes remain in the separate ROKS CSV.
The second prints full and reduced QUEST1 MAEs and writes
`output/tables/quest1_outlier_summary.csv` and `quest1_outlier_errors.csv`.
Ethylene, Cyclopropene and Streptocyanine-C1 are excluded identically from every
comparison. This is descriptive sensitivity analysis, not statistical rejection;
the full benchmark remains primary.

## Data and regeneration

To rebuild the benchmark CP2K excitation CSV from the included outputs:

```bash
uv run python scripts/parse_cp2k.py
```

See the [ROKS data guide](data/roks-st-gap/README.md) for rebuilding
`data/csv/roks_st_gap.csv` and regenerating frozen-orbital inputs.
[Reference sources](data/references/README.md) and
[PySCF provenance](data/pyscf/README.md) are documented separately.

Wavefunctions (`*.wfn`) are excluded from Git and retained for the Zenodo
archive. They are needed to regenerate frozen-orbital inputs or repeat density
checks, but not to reproduce figures and tables from the supplied CSV files.
