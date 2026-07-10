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

`magmoms` are the site-projected magnetic moments of the two atoms (`LORBIT=11`, in μB) and `spin_candidates` records which spin-ladder rung won the per-distance energy minimum (an even `NUPDOWN` integer as string, or `afm` for the broken-symmetry candidate). Both serve two purposes: discontinuities in atom-wise moments vs distance flag SCF spin-state hops that the total (pinned) moment cannot show, and spin-aware MLIPs (which take magnetic moments or spin as input) can be compared against the reference spin state directly. A few points (6 of 9113) have `null` magmoms where the OUTCAR block was unavailable.

## DFT Setup

Each candidate curve was produced with VASP 6 using Materials Project (MP24)-style static input settings and real-space PAW potentials for PBE and r2SCAN.

The PBE 6.4 PAW labels used for `H...U` were the `MP24StaticSet` labels:

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

Fixed-oxidation-state lanthanide potentials such as `Ho_3` freeze f-electrons in the core with the occupancy of the trivalent ion, so they are not appropriate for neutral Ho₂. We instead use the `_h` potentials (`Pr_h` through `Yb_h`), which include f-electrons in the valence. This choice also makes SCF convergence difficult for heavy lanthanides (see caveats below). The current Lu curve uses `Lu_3`, which freezes 4f¹⁴ but retains the neutral atom's 5d¹6s² bonding electrons (plus 5p⁶ semicore). Matching the neutral f-shell occupancy does not by itself establish transferability to Lu₂, so we are validating it against the standard `Lu` potential with 4f¹⁴ in the valence. For Rh, we use `Rh_pv`, which includes semicore p states in the valence and does not impose an oxidation state.

Important per-point settings:

- `ISMEAR=0`, `SIGMA=0.05`.
- `KSPACING=0.6` (Gamma-only in this large cell). A 2x2x2 $k$-grid spot check on Na/PBE and Al/PBE changed final energies by 2.3 meV mean absolute over 40 matched points; Na's diffuse 6 Å tail (the largest separation) reached 23 meV, while Al stayed below 6.2 meV.
- `ISPIN=2` with fixed `NUPDOWN` for every spin candidate and explicit `MAGMOM` initialization where needed.
- `NELM=250` for difficult stretched or repulsive geometries.
- `LWAVE=.TRUE.` so neighboring geometries can reuse the wavefunction.
- `LCHARG`, `LAECHG`, `LVTOT`, and `LELF` disabled to avoid roughly 180 MB of unused grid output per point; `LORBIT=11` to record site-projected magnetic moments for every point (the spin-hop diagnostic and spin-aware-MLIP reference described above).
- With a seed wavefunction: `ISTART=1`, `ICHARG=0`; otherwise cold start with `ISTART=0`, `ICHARG=2`.

The distance grid is 50 geometrically spaced separations from `0.8 * r_cov` to 6 Å. Going deeper into the repulsive wall caused severe PAW augmentation-sphere overlap and hundreds-of-eV energies, which made the reference less useful for MLIP scoring.

The 15 Å box was chosen after checking the vacuum error/cost tradeoff. A 13 Å box left roughly 0.15 eV periodic-image errors for diffuse large-distance atoms; 15 Å reduced this to roughly 0.03 eV; 18 Å was closer to converged but about 2.7x more expensive.

## Spin Ladder

The reference is not a single fixed-spin sweep. SCF branch hops between spin multiplicities can create cliffs and spikes, so each element uses a spin ladder:

- even `NUPDOWN` values `0, 2, ...` up to the larger of the known molecular ground-state spin and `2 *` the free-atom Hund moment;
- one extra rung for elements whose dimer spin state is uncertain;
- an additional broken-symmetry AFM candidate (`NUPDOWN=0`, `MAGMOM=+m -m`) for every open-shell atom.

The candidate map is a JSON object `{ element: candidate[] }` where each candidate is an even integer `NUPDOWN` value or the string `afm`. It defines which spin candidates are merged for each element and functional. Raw candidate curve directories follow this naming convention:

```text
<El>_<xc>_n<NUPDOWN>/curve.json
<El>_<xc>_n0afm/curve.json
```

## Warm Starts

For most fixed-spin candidates, the sweep cold-started near the bonded anchor (nearest `2 * r_cov`) and then warm-started in both directions. Inward points seeded from the next larger separation; outward points seeded from the next smaller separation. The constant box keeps FFT grids compatible, so `WAVECAR` or `vaspwave.h5` can be reused.

AFM candidates are different: they anchor at the largest separation and sweep inward. This preserves the broken-symmetry free-atom solution, which otherwise tends to collapse if the chain is started near the bonded geometry.

After the main sweep, a conservative reconciliation pass retried obvious warm-start branch traps. Interior energy bumps of at least 0.3 eV were retried from the lower-energy neighbor, with hard caps on the number of passes and retries to avoid f-electron thrashing.

## Building The Bundled Reference

The checked-in builder is `scripts/evals/build_diatomic_reference.py`.

For each element and functional, it:

1. reads the per-element spin state candidate map;
1. loads finite-energy, finite-force points from the spin-candidate curves;
1. merges candidate curves by taking the per-distance minimum energy;
1. applies the checked-in postprocess in `matbench_discovery.metrics.diatomics.reference` to remove isolated SCF artifacts (see [Postprocessing](#postprocessing));
1. checks the merged PECs for branch jumps and missing candidate coverage;
1. writes the compact, public, version-controlled site artifact to `site/src/lib/diatomics-dft.json.gz`;
1. rebuilds model diatomic metrics and confirms the DFT curves display on the Diatomics page.

The artifact contains 92 PBE and 92 r2SCAN homonuclear entries. Every raw spin-candidate curve now meets the 45-point completeness threshold. Merged curves can be shorter because postprocessing drops invalid points.

The quality checks are diagnostic and do not alter endpoints. `count_dissociation_tail_jumps` inspects adjacent energy steps among up to the final three merged points, flags steps of at least 0.1 eV, and ignores the short-range repulsive wall. `reference-quality.json` records the number of flagged steps as `tail_jumps`; the `tail_jump_pairs` summary counts each affected element-functional pair once. The current artifact flags 7 PBE and 12 r2SCAN pairs. The clearest known endpoint spin-branch hops are `O/r2SCAN` and `Rh/r2SCAN` near 6 Å. `Ho/r2SCAN`, `Er/r2SCAN`, and similar heavy lanthanides and actinides remain branch-trapped and jumpy despite all postprocessing and checks. The heavy-element failures seem irrecoverable, at least with current VASP 6.4 pseudo-potentials.

## Postprocessing

The postprocess is deliberately narrow. It is meant to fix isolated SCF artifacts without smoothing away real physics like spin crossovers. Each filter targets a failure mode no smooth adiabatic PEC can produce; anything less clear-cut is published as computed (per-point magmoms and spin candidates make residual branch flicker visible rather than hidden).

It handles three cases:

- variationally collapsed SCF points, dropped per spin candidate before merging: mostly f-electron heavy lanthanides where single points fall tens to thousands of eV below their own branch (no physical PEC well is more than a few eV deep, and such points would otherwise always win the per-distance minimum);
- isolated severe (&ge;3 eV) spin-branch drops in the merged curve, replaced by the smoother neighboring branch (a real spin crossover stays lower over a range of separations; a single-point notch is an SCF failure);
- isolated upward energy bumps (&ge;0.1 eV single-grid-point local maxima, which no smooth PEC has), dropped.

An earlier revision also substituted short, near-degenerate (&le;0.2 eV) "spin-branch islands" with the surrounding branch for cosmetic smoothness. That step was removed: it edited far more points than all other filters combined (427 of 605 edits) for negligible smoothness gain, and replacing the computed adiabatic minimum with a higher-energy branch is hand-curation, not physics.

Tests for this logic live in `tests/metrics/diatomics/test_diatomics_reference.py` and `tests/metrics/diatomics/test_build_diatomic_reference.py`.

## Acknowledgements

Thanks to [Andrew S. Rosen](https://cbe.princeton.edu/people/andrew-rosen) [[Google Scholar](https://scholar.google.com/citations?user=lHBjgLsAAAAJ&hl=en)] for guidance on how to obtain high quality dimer energy curves!
