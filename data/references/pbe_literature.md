# Published PBE calculations

`pbe_literature.csv` contains published orbital-optimized PBE excitation
energies, **not TBEs**. These values are kept separate from `excitations.csv`.
The literature data are transcribed to the published precision (two decimals,
nominal rounding of ±0.005 eV). Dataset/molecule identifiers match the CP2K CSV.

## Sources and target selection

- **QUEST1:** Ivanov, Levi, Jónsson and Jónsson, *Method for Calculating Excited
  Electronic States Using Density Functionals and Direct Orbital Optimization
  with Real Space Grid or Plane-Wave Basis Set*, JCTC 17, 5034–5049 (2021).
  https://doi.org/10.1021/acs.jctc.1c00157 . PBE triplets are from Table 1,
  printed page 5042; spin-purified PBE singlets from Table 2, page 5043.
  Seven matched systems: acetaldehyde, ammonia, diazomethane, ethylene,
  methanimine, thioformaldehyde, and water (14 excitation energies).
- **QUEST6:** Selenius, Sigurdarson, Schmerwitz and Levi, *Orbital-Optimized Versus
  Time-Dependent Density Functional Calculations of Intramolecular Charge
  Transfer Excited States*, JCTC 20, 3809–3822 (2024).
  https://doi.org/10.1021/acs.jctc.3c01319 . All eleven selected singlets are
  from Table 1, page 3814, the **orbital optimized / PBE / ΔE** column, not TDDFT.

For QUEST1, the paper's formaldehyde and ketene entries correspond to different
excited states from our selected targets, so they are not compared. Cyclopropene,
nitrosomethane and streptocyanine-C1 are absent from those paper tables. The
ammonia match follows excitation character and the reference-energy pair;
published symmetry labels depend on the coordinate convention. This comparison
does not independently validate that convention or spectroscopic assignment.

## Numerical conventions and limits

The comparison uses the current CP2K UKS/PBE/aug-cc-pVTZ CSV. UKS singlet energies
are spin purified as `2 E(mixed) - E(triplet)`. Signed differences are CP2K minus
the published PBE values. MAD and MaxAD describe differences between calculations,
not errors against TBEs. QUEST1 summaries are given separately for singlets and
triplets and jointly over all 14 energies; QUEST6 covers 11 singlets.

The literature protocols differ from our all-electron GAPW/aug-cc-pVTZ protocol:
Ivanov uses GPAW frozen-core PAW with a finite-difference real-space grid
(0.2 Å spacing) and direct optimization with MOM. Selenius uses GPAW frozen-core
PAW, primitive aug-cc-pVDZ Gaussians plus single-zeta numerical atomic orbitals,
a 0.15 Å grid, and freeze-and-release direct optimization. Geometry sources
are the corresponding QUEST benchmarks, but no independent atom-by-atom
literature geometry comparison was performed in the reviewed extraction.
Representation, basis, grid, cell and optimization effects are not isolated.
Energy agreement supports consistency of the nominal targets; it does not prove
state assignment.

The published values and table/page locations were transferred from the reviewed
literature extraction on 2026-11-13. Historical CP2K values were not copied:
`table_pbe_literature.py` always joins to the current submission CSV.
