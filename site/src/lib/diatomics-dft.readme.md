# Homonuclear Diatomic DFT Reference Curves

`diatomics-dft.json.gz` contains VASP reference potential-energy curves for homonuclear dimers `X2` with `X = H...U` (`Z = 1...92`) at PBE and r2SCAN. The site uses this file for the DFT overlays on the Diatomics page and the Python metrics code uses the PBE subset for reference-relative model metrics.

After decompression, the schema is:

```ts
{
  [functional: string]: {
    [formula: string]: {
      distances: number[]
      energies: number[]
      forces: number[][][] // [point][atom][x|y|z]
      magmoms: ([number, number] | null)[] // per-point site-projected moments (μB)
      spin_candidates: (string | null)[] // winning spin candidate, "<NUPDOWN>" or "afm"
    }
  }
}
```

`functional` is `PBE` or `r2SCAN`; `formula` is e.g. `O-O`. `distances`, `energies`, `forces`, `magmoms`, and `spin_candidates` share the same first axis over curve points. Energies are final sigma-to-0 VASP energies in eV; forces are per-point, per-atom Cartesian 3-vectors in eV/Å with axis order `[point][atom][x|y|z]`.

`magmoms` are the `LORBIT=11` site projections for the two atoms (in μB), and `spin_candidates` records which constrained-magnetization candidate won the per-distance energy minimum (an even `NUPDOWN` integer as a string, or `afm` for the broken-symmetry candidate). The site projections are qualitative and depend on the PAW projectors; they exclude interstitial magnetization and are not unique atomic spin labels. Abrupt changes can nevertheless help diagnose changes in the converged electronic branch that the fixed total `NUPDOWN` cannot show. These values can aid analysis of spin-aware MLIPs but are directly comparable only to models using the same VASP projection convention. A few points (6 of 9113) have `null` magmoms where the OUTCAR block was unavailable.

## DFT Setup

Each candidate curve was produced with VASP 6 using Materials Project (MP24)-style static input settings. Both PBE and r2SCAN calculations used PBE_64 PAW pseudo-potentials.

The resulting PBE_64 PAW labels for `H...U`, including the overrides discussed below, were:

```yaml
H: H
He: He
Li: Li_sv
Be: Be_sv
B: B
C: C
N: N
O: O
F: F
Ne: Ne
Na: Na_pv
Mg: Mg_pv
Al: Al
Si: Si
P: P
S: S
Cl: Cl
Ar: Ar
K: K_sv
Ca: Ca_sv
Sc: Sc_sv
Ti: Ti_pv
V: V_pv
Cr: Cr_pv
Mn: Mn_pv
Fe: Fe_pv
Co: Co
Ni: Ni_pv
Cu: Cu_pv
Zn: Zn
Ga: Ga_d
Ge: Ge_d
As: As
Se: Se
Br: Br
Kr: Kr
Rb: Rb_sv
Sr: Sr_sv
Y: Y_sv
Zr: Zr_sv
Nb: Nb_pv
Mo: Mo_pv
Tc: Tc_pv
Ru: Ru_pv
Rh: Rh_pv
Pd: Pd
Ag: Ag
Cd: Cd
In: In_d
Sn: Sn_d
Sb: Sb
Te: Te
I: I
Xe: Xe_GW
Cs: Cs_sv
Ba: Ba_sv_GW
La: La
Ce: Ce
Pr: Pr_h
Nd: Nd_h
Pm: Pm_h
Sm: Sm_h
Eu: Eu
Gd: Gd
Tb: Tb_h
Dy: Dy_h
Ho: Ho_h
Er: Er_h
Tm: Tm_h
Yb: Yb_h
Lu: Lu_3
Hf: Hf_pv
Ta: Ta_pv
W: W_sv
Re: Re_pv
Os: Os_pv
Ir: Ir
Pt: Pt
Au: Au
Hg: Hg
Tl: Tl_d
Pb: Pb_d
Bi: Bi
Po: Po_d
At: At
Rn: Rn
Fr: Fr_sv
Ra: Ra_sv
Ac: Ac
Th: Th
Pa: Pa
U: U
```

Lanthanide `_3` PAW datasets use a valence configuration intended for trivalent environments and freeze the f-electrons in the core. They do not strictly fix an oxidation state, but that frozen-core approximation likely transfers poorly to a neutral dimer. We therefore use datasets with f-electrons in the valence: `_h` variants for Pr-Sm and Tb-Yb, and the standard `Eu` and `Gd` datasets. Their larger valence spaces cost more, and the heavy-lanthanide SCF branches proved difficult to converge with the present protocol (see caveats below). The current Lu curve still uses `Lu_3` (`ZVAL=9`). Although its frozen 4f¹⁴ shell has the same nominal occupancy as neutral Lu, that alone does not establish transferability to Lu₂ (analysis for `Lu` PSP with `ZVAL=25` which treats 4f¹⁴ explicitly ongoing). For Rh, `Rh_pv` includes the 4p semicore states in the valence.

Important per-point settings:

- `ISMEAR=0`, `SIGMA=0.05`.
- `KSPACING=0.6` (Gamma-only in this large cell). A 2x2x2 $k$-grid spot check on Na/PBE and Al/PBE changed final energies by 2.3 meV mean absolute over 40 matched points; Na's diffuse 6 Å tail (the largest separation) reached 23 meV, while Al stayed below 6.2 meV.
- `ISPIN=2` with fixed `NUPDOWN` for every spin candidate and explicit `MAGMOM` initialization where needed.
- `NELM=250` for difficult stretched or repulsive geometries.
- `LWAVE=.TRUE.` so neighboring geometries can reuse the wavefunction.
- `LCHARG`, `LAECHG`, `LVTOT`, and `LELF` disabled to avoid roughly 180 MB of unused grid output per point; `LORBIT=11` to request the site-projected magnetic moments used as the branch diagnostic described above.
- With a seed wavefunction: `ISTART=1`, `ICHARG=0`; otherwise cold start with `ISTART=0`, `ICHARG=2`.

The distance grid is 50 geometrically spaced separations from `0.8 * r_cov` to 6 Å. The repulsive-wall metric uses the full DFT range and compares wall radii at 1, 5, 10, 20, 50, and 100 eV above the well wherever the reference reaches those energies. A prediction that does not reach a supported threshold receives the full reference-radius error instead of having that threshold omitted. General energy, force, and smoothness metrics start at `0.9 * r_cov` so the hundred-eV endpoint does not dominate every metric. We did not run DFT below `0.8 * r_cov` because PAW augmentation-sphere overlap increases at shorter distances. No reference data exist below that cutoff, so those separations are not scored.

The 15 Å box was chosen after checking the vacuum error/cost tradeoff. For the tested diffuse, large-separation atoms, such as Na₂ at 6 Å, energies at 13 and 15 Å differed from the corresponding 18 Å results by roughly 0.15 and 0.03 eV, respectively; 18 Å cost about 2.7x more than 15 Å. These differences estimate residual periodic-image effects but are not an extrapolation to infinite vacuum.

## Spin Ladder

The reference is not a single fixed-`NUPDOWN` sweep. Different constrained collinear magnetizations can be lowest at different separations, and the SCF can also converge to different orbital solutions within one sector. Each element therefore uses a candidate ladder:

- even `NUPDOWN` values `0, 2, ...` up to the larger of the reported molecular `2S` value and `2 *` the free-atom Hund moment;
- one extra rung for elements whose dimer spin state is uncertain;
- an additional broken-symmetry AFM candidate (`NUPDOWN=0`, `MAGMOM=+m -m`) for every open-shell atom.

The candidate map is a JSON object `{ element: candidate[] }` where each candidate is an even integer `NUPDOWN` value or the string `afm`. It defines which spin candidates are merged for each element and functional. Raw candidate curve directories follow this naming convention:

```text
<El>_<xc>_n<NUPDOWN>/curve.json
<El>_<xc>_n0afm/curve.json
```

## Warm Starts

For most fixed-`NUPDOWN` candidates, the sweep cold-started near the bonded anchor (nearest `2 * r_cov`) and then warm-started in both directions. Inward points seeded from the next larger separation; outward points seeded from the next smaller separation. The constant box keeps FFT grids compatible, so `WAVECAR` or `vaspwave.h5` can be reused.

AFM candidates are different: they anchor at the largest separation and sweep inward. This is intended to retain a broken-symmetry, free-atom-like solution and reduced the observed tendency to collapse when the chain was started near the bonded geometry; it does not guarantee that the same electronic branch is retained throughout.

After the main sweep, a conservative reconciliation pass retried suspected warm-start branch traps. Interior energy bumps of at least 0.3 eV were treated as possible artifacts and retried from the lower-energy neighbor, with hard caps on the number of passes and retries to avoid f-electron thrashing.

## Building The Bundled Reference

The checked-in builder is `scripts/evals/build_diatomic_reference.py`.

For each element and functional, it:

1. reads the per-element constrained-magnetization candidate map;
1. loads finite-energy, finite-force points from the spin-candidate curves;
1. merges candidate curves by taking the per-distance minimum energy;
1. applies the checked-in postprocess in `matbench_discovery.metrics.diatomics.reference` to filter points classified by local heuristics as likely SCF artifacts (see [Postprocessing](#postprocessing));
1. checks the merged PECs for energy discontinuities and missing candidate coverage;
1. writes the compact, public, version-controlled site artifact to `site/src/lib/diatomics-dft.json.gz`;
1. rebuilds model diatomic metrics and confirms the DFT curves display on the Diatomics page.

The artifact contains 92 PBE and 92 r2SCAN homonuclear entries. Every raw spin-candidate curve now meets the 45-point completeness threshold. Merged curves can be shorter because postprocessing drops invalid points.

The quality checks are diagnostic and do not alter endpoints. `count_dissociation_tail_jumps` inspects adjacent energy steps among up to the final three merged points, flags steps of at least 0.1 eV, and ignores the short-range repulsive wall. `reference-quality.json` records the number of flagged steps as `tail_jumps`; the `tail_jump_pairs` summary counts each affected element-functional pair once. The current artifact flags 7 PBE and 12 r2SCAN pairs. The largest known endpoint discontinuities associated with changes in the selected candidate or projected moments include `O/r2SCAN` and `Rh/r2SCAN` near 6 Å. `Ho/r2SCAN`, `Er/r2SCAN`, and several other heavy lanthanide and actinide curves still show large discontinuities or strong branch sensitivity. The present spin ladder and warm-start retries did not resolve those cases; this does not show that alternative initializations, PAW datasets, relativistic treatments, or electronic-structure methods could not.

## Postprocessing

The postprocess is deliberately narrow. It applies local heuristics to points that are strongly suggestive of SCF artifacts while trying not to remove physical features such as barriers or state crossings. These classifications are not proofs: a sufficiently narrow physical feature could resemble a one-grid-point artifact. Anything that does not meet the explicit thresholds is published as computed, and the per-point projected moments and selected candidates expose residual branch sensitivity.

It handles three cases:

- candidate points that fall tens to thousands of eV below their neighboring branch, seen mostly for f-electron heavy lanthanides, are dropped as numerically implausible discontinuities before merging because they would otherwise always win the per-distance minimum;
- isolated severe (&ge;3 eV) drops in the merged curve are replaced by the neighboring candidate when the lower state appears at only one sampled separation; this treats a one-point notch as a likely branch failure
- isolated upward bumps (&ge;0.1 eV single-grid-point local maxima) are dropped as likely SCF artifacts, not because physical PECs are forbidden from having barriers or multiple wells.

An earlier revision also substituted short, near-degenerate (&le;0.2 eV) "spin-branch islands" with the surrounding branch for cosmetic smoothness. That step was removed: it edited far more points than all other filters combined (427 of 605 edits) for negligible smoothness gain and imposed a stronger hand-curation assumption by replacing the computed minimum with a higher-energy branch.

Tests for this logic live in `tests/metrics/diatomics/test_diatomics_reference.py` and `tests/metrics/diatomics/test_build_diatomic_reference.py`.

## Acknowledgements

Thanks to [Andrew S. Rosen](https://cbe.princeton.edu/people/andrew-rosen) [[Google Scholar](https://scholar.google.com/citations?user=lHBjgLsAAAAJ&hl=en)] for guidance on how to obtain high quality dimer energy curves!
