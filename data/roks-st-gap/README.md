# Frozen ROKS/ST-gap data

The original PBE–Skala study and the new r2SCAN/PBE0 evaluations are separate
snapshots. All use the eight selected ST-gap molecules, fixed nuclear geometries,
and aug-cc-pVTZ. No ground-state energies enter this analysis.

For evaluating functional F on donor G orbitals, M and TS use the same
ROKS-singlet orbitals; TT uses the separately optimized restricted-triplet
orbitals. The reconstructed singlet is 2M−TS; the common-orbital splitting is
2(M−TS); the triplet correction is TS−TT; the reported gap is 2M−TS−TT.
Cross-evaluated TT is frozen on G's triplet orbitals, not reoptimized with F.
The PBE–Skala functional/orbital partition is path-dependent and distinct from
the same-functional triplet-relaxation test.

## Contents

- `original-pbe-skala/`: two disjoint November 6 campaigns, 96 probes, 32
  functional/orbital combinations, preserving the figure's original endpoints.
- `r2scan-pbe0/`: November 13 same-functional campaign, with donor reproduction
  checks using the separately optimized triplet wavefunctions.
- Each campaign has `sources/<functional>/<molecule>/`, `trials/<id>/`, a
  manifest and execution records. Historical absolute paths are provenance
  strings only; the parser resolves included files relative to each campaign.
- Each probe includes its input/output, pre-update spin density dump,
  validation record, runtime information and seed WFN. Sources retain donor
  singlet/triplet inputs, outputs and WFNs. Donor inputs record historical
  restart chains; frozen regeneration starts from the included final WFNs.
- Wavefunctions (`*.wfn` and `*.wfn.bak*`) and density dumps (`density.Log`
  and `density.Log.gz`) are ignored by Git but retained for the Zenodo archive.
  The parser can verify stored validation records without WFNs. Repeating
  numerical density checks or generating inputs requires the WFN archive.

The original endpoints differ from the currently adopted benchmark gaps by at
most 0.000020 eV; `original-pbe-skala/endpoint_comparison.csv` records every
comparison. They remain separate snapshots. The new parser retains full printed
energy precision, so differences from the earlier iteration-row CSV are only a
few nano-eV and do not affect reported rounding or inversion counts.

`cp2k-data/` contains the basis/potential files and their SHA256 hashes. These
match the original and current runtime records. When reproducing an input,
make these files available through `CP2K_DATA_DIR` or the working directory.
Skala evaluations require the compatible CP2K/Skala installation and the
Skala 1.1 rev1 model identified by the runtime hash; model/software binaries are
not bundled here. Archived launchers record the original site configuration.

## Commands (repository root)

Generate new inputs without changing an existing campaign:

```bash
uv run python scripts/generate_roks_st_gap.py \
  --source data/roks-st-gap/r2scan-pbe0/roks-st-gap-2026-11-13/sources \
  --destination /tmp/roks-new-inputs \
  --pairs r2scan@r2scan pbe0@pbe0
```

`--pairs pbe@pbe skala@pbe pbe@skala skala@skala` prepares the original comparison
matrix when supplied its original donor snapshots. Do not use current production
PBE/Skala donors to reconstruct the historical figure. The archived original
campaigns split the molecules into two groups: the `1224UTC` campaign contains
Heptazine, Cyclazine, Molecule6 and Molecule8; `1317UTC` contains Molecules3–5
and Molecule7. Use each campaign's `sources/` as `--source`, the four pairs above,
and the matching `--molecules` list to regenerate these snapshots separately.
Generator outputs are not automatically submitted. The archived launch script
records the actual RoSi environment and resources, rather than prescribing a
portable cluster configuration. Its original 16-rank/64GB pilot allocation
failed for PBE0; the successful PBE0 allocations and exact array submissions
are recorded in `r2scan-pbe0/roks-st-gap-2026-11-13/submission.json`.

Parse selected complete campaigns:

```bash
uv run python scripts/parse_roks_st_gap.py \
  --campaign data/roks-st-gap/original-pbe-skala/roks-skala-frozen-2026-11-06-1224UTC \
  --campaign data/roks-st-gap/original-pbe-skala/roks-skala-frozen-2026-11-06-1317UTC \
  --campaign data/roks-st-gap/r2scan-pbe0/roks-st-gap-2026-11-13 \
  --output data/csv/roks_st_gap.csv
```

Add `--validate-densities` to repeat seed-density validation using the WFN files
and density dumps from the Zenodo archive.
The full inventory, job outcomes, numerical checks and remaining qualifications
are in [VALIDATION.md](VALIDATION.md).

The scripts require only Python and NumPy and do not import the development
repository. Full printed total energies are retained; 27.211396132 eV/Ha is used.
The adjacent validation JSON reports means, inversion counts, sign changes,
donor-energy reproduction and non-XC residuals including GAPW contributions.
PBE0 XC includes exact exchange.

These are deliberately stopped one-step UKS determinant evaluations of ROKS
orbitals. Their first pre-update energies are accepted only with separate
frozen-density validation; they are not converged UKS benchmarks. Numerical
reproduction and orbital correspondence do not establish spectroscopic state
assignment. In particular, a near-zero gap sign is not a grid/basis convergence
claim. The archive preserves the existing study's qualifications.
