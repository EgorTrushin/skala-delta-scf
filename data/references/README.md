# Reference excitation energies

`excitations.csv` contains the theoretical best estimates (TBEs) adopted in the
manuscript: 12 QUEST1 systems, eight ST-gap systems, and 11 QUEST6 systems.
All energies are in eV. Join to calculated results on `(dataset, molecule)`;
identifiers match the calculation directories and `cp2k_excitations.csv`.

| Dataset | Source and adopted values |
| --- | --- |
| QUEST1 | Loos et al., *A Mountaineering Strategy to Excited States: Highly Accurate Reference Energies and Benchmarks*, JCTC 14, 4360–4379 (2018), corrected TBEs in Table 6. DOI: https://doi.org/10.1021/acs.jctc.8b00406 |
| ST_gap | Vertical S1/T1 TBEs in Supporting Information Table S1, p. S2. DOI: https://doi.org/10.1021/acs.jpclett.3c03042 |
| QUEST6 | Loos et al., *Reference Energies for Intramolecular Charge-Transfer Excitations*, JCTC 17, 3666–3686 (2021), final aug-cc-pVQZ TBEs in Table 1. DOI: https://doi.org/10.1021/acs.jctc.1c00226 |

The QUEST1 ethylene singlet TBE is **7.92 eV** (not the earlier 7.89 eV value).
The selected ethylene and cyclopropene singlets are not their lowest singlet
excited states. Reference gaps are calculated as `singlet_eV - triplet_eV`
for the selected states; no separate gap column is stored.

QUEST6 has no reference triplet energies for this benchmark; empty triplet
fields mean unavailable, not zero. Its adopted QZ TBEs are used to assess the
TZ calculations. The same reference values are used for all functionals,
methods and basis sets. These published targets do not establish assignment
of the calculated electronic states.

The values were transferred from the reviewed manuscript reference inventory
on 2026-11-13 and checked against the adopted reference data: 31 systems and
51 excitation TBEs. Every system in the initial 348-row CP2K export has a
matching reference. The reference CSV is a supplied input, not regenerated
from calculation outputs.

## Attribution

The QUEST data retain their Creative Commons Attribution-ShareAlike 4.0
International (CC BY-SA 4.0) terms, separately from the analysis code:
https://creativecommons.org/licenses/by-sa/4.0/ .
QUESTDB: https://github.com/pfloos/QUESTDB ; dataset snapshot:
https://doi.org/10.5281/zenodo.15671384 . The source publication for each
reference row is also recorded in the CSV. This file contains a selected
subset with directory-compatible molecule identifiers and the adopted
ethylene value noted above.

## State metadata

`states.json` records the published state symmetries, excitation character,
QUEST6 reference root numbers, ST-gap compound numbering, molecular charge,
and source-table locations used by the first three SM tables. Entries retain
the reviewed table order and map to the energy CSV through dataset and molecule
identifiers. Reference labels do not establish calculated-state assignments.
For QUEST1 the table character column describes the singlet; the metadata
retain the distinct published triplet character where applicable (ketene).
QUEST6 root numbers refer to the publication's LR-CCSD/cc-pVTZ ordering.
The metadata were transferred from the reviewed reference inventory on
2026-11-13; all three generated tables match the current SM source.

## Published PBE results (separate from TBEs)

`pbe_literature.csv` contains 25 published PBE excitation energies used for
the implementation/target consistency comparison. See
[pbe_literature.md](pbe_literature.md) for citations, selected states and
protocol qualifications. These values are not used to calculate the
main-text MAEs against theoretical best estimates.
