# Homonuclear Diatomic DFT Reference Curves

`diatomics-dft.json.gz` contains reference potential-energy curves for homonuclear dimers `X2` with `X = H...U` (`Z = 1...92`) at PBE and r2SCAN. The site uses this file for the DFT overlays on the Diatomics page and the Python metrics code uses the PBE subset for reference-relative model metrics.

Download the curves from [Figshare](https://figshare.com/files/68541277), or use `DataFiles.diatomics_dft_reference.path` to download, verify, and cache them in Python. The bundled copy keeps site builds independent of Figshare availability. When regenerating it, publish the updated file and update the registry URL and checksum together.

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

`functional` is `PBE` or `r2SCAN`; `formula` is e.g. `O-O`. `distances`, `energies`, `forces`, `magmoms`, and `spin_candidates` share the same first axis over curve points. Energies are final sigma-to-0 energies in eV; forces are per-point, per-atom Cartesian 3-vectors in eV/Å with axis order `[point][atom][x|y|z]`.

`magmoms` are the `LORBIT=11` site projections for the two atoms (in μB), and `spin_candidates` records which constrained-magnetization candidate won the per-distance energy minimum (an even `NUPDOWN` integer as a string, or `afm` for the broken-symmetry candidate). The site projections are qualitative and depend on the PAW projectors; they exclude interstitial magnetization and are not unique atomic spin labels. Abrupt changes can nevertheless help diagnose changes in the converged electronic branch that the fixed total `NUPDOWN` cannot show. These values can aid analysis of spin-aware MLIPs but are directly comparable only to models using the same projection convention.

## DFT Setup

Each candidate curve was computed using Materials Project (MP24)-style static settings with one deliberate change: `ISYM=0` (see [Symmetry Breaking](#symmetry-breaking-isym0)). Both PBE and r2SCAN calculations used PBE_64 PAW pseudo-potentials.

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
Lu: Lu
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

Lanthanide `_3` PAW datasets use a valence configuration intended for trivalent environments and freeze the f-electrons in the core. They do not strictly fix an oxidation state, but that frozen-core approximation seemed dangerous as it might transfer poorly to a neutral dimer (tl;dr empirically seems fine, saw little between e.g. `Lu_3` and `Lu`). To be safe, we nonetheless opted for datasets with f-electrons in the valence: `_h` variants for Pr-Sm and Tb-Yb, and the standard `Eu`, `Gd`, and `Lu` datasets. The standard `Lu` potential (`ZVAL=25`) treats its 4f¹⁴ shell explicitly and replaces the earlier `Lu_3` (`ZVAL=9`) curves (though as said, the difference was minor). These larger valences cost more, and the heavy-lanthanide SCF branches need the occupation seeding described under [Occupation Seeding](#occupation-seeding) to reach their ground states. For Rh, `Rh_pv` includes the 4p semicore states in the valence.

Important per-point settings:

- `ISYM=0`: no density symmetrization.
- `ISMEAR=0`, `SIGMA=0.05`.
- `KSPACING=0.6` (Gamma-only in this large cell). A 2x2x2 $k$-grid spot check on Na/PBE and Al/PBE changed final energies by 2.3 meV mean absolute over 40 matched points; Na's diffuse 6 Å tail (the largest separation) reached 23 meV, while Al stayed below 6.2 meV.
- `ISPIN=2` with fixed `NUPDOWN` for every spin candidate and explicit `MAGMOM` initialization where needed.
- `NELM=250` for difficult stretched or repulsive geometries.
- `LWAVE=.TRUE.` so neighboring geometries can reuse the wavefunction.
- `LCHARG`, `LAECHG`, `LVTOT`, and `LELF` disabled to avoid roughly 180 MB of unused grid output per point; `LORBIT=11` to request the site-projected magnetic moments used as the branch diagnostic described above.
- Warm points (a neighboring geometry's wavefunction as the start): `ISTART=1`, `ICHARG=0`, `ALGO=Normal` (Davidson).
- Cold points (the first point of a chain, or a retry): `ISTART=0`, `ICHARG=2`, `ALGO=All` with `NELMDL=-8`. Open-shell atoms and dimers stall in density mixing between near-degenerate levels, which the all-band minimizer avoids; the 8 non-self-consistent delay steps also hold an occupation seed (see [Occupation Seeding](#occupation-seeding)) before the SCF relaxes freely.
- An SCF that stalls is restarted from scratch with the all-band minimizer, and one whose warm start aborts without an energy (the PAW one-center check at the shortest separations) is retried cold. Each point gets at most 1,000 SCF steps over all attempts; a point that reaches the cap is dropped and its chain continues from the wavefunction that point started from. Points that still do not converge are dropped.

The distance grid is 50 geometrically spaced separations from `0.8 * r_cov` to 6 Å. The repulsive-wall metric uses the full DFT range and compares wall radii at 1, 5, 10, 20, 50, and 100 eV above the well wherever the reference reaches those energies. A prediction that does not reach a supported threshold receives the full reference-radius error instead of having that threshold omitted. General energy, force, and smoothness metrics start at `0.9 * r_cov` so the hundred-eV endpoint does not dominate every metric. We did not run DFT below `0.8 * r_cov` because PAW augmentation-sphere overlap increases at shorter distances. No reference data exist below that cutoff, so those separations are not scored.

The 15 Å box was chosen after checking the vacuum error/cost tradeoff. For the tested diffuse, large-separation atoms, such as Na₂ at 6 Å, energies at 13 and 15 Å differed from the corresponding 18 Å results by roughly 0.15 and 0.03 eV, respectively; 18 Å cost about 2.7x more than 15 Å. These differences estimate residual periodic-image effects but are not an extrapolation to infinite vacuum.

## Symmetry Breaking (`ISYM=0`)

The first version of this reference (Figshare file [68541277](https://figshare.com/files/68541277)) used the MP24 default `ISYM=2`, which symmetrizes the charge density under the point group of the dimer in its box. That excluded the orbitally polarized ground states of open-shell dimers and left the heavy-lanthanide curves with tens-of-eV discontinuities. The current version disables symmetrization and is compared against the first version below, for these reasons:

- **Symmetrization forbids the ground state.** For open d, f and p shells the lowest-energy solution is often orbitally polarized, i.e. it breaks the dimer's `D∞h` symmetry. `ISYM=2` excludes such states by construction, so the symmetric curves describe a constrained excited state whose energy depends on how the code handles symmetry rather than on the physics. Disabling it lowers free f atoms by 63-110 meV (PBE Tm, Er, Ho), d-block dimers by 20-400 meV, and r2SCAN open-shell dimers by up to ~1 eV (e.g. O₂ by 0.8 eV). Closed-shell dimers are unchanged to within 0.1 meV.
- **A variational reference should be the lowest state found.** Every point is the minimum over all spin candidates, seeds and sweep directions, and lies at or below the first version's `ISYM=2` value at nearly every shared distance (PBE: 21 pairs agree within 5 meV everywhere, 54 lower everywhere, 16 mixed, 1 marginally higher by at most 34 meV; r2SCAN: 4 agree, 69 lower, 19 mixed, 0 higher).
- **Symmetric SCF produced broken curves.** The `ISYM=2` Ho, Er and Tm curves had tens-of-eV discontinuities and Dy/Pr eV-scale kinks. With `ISYM=0` plus occupation seeding all 184 curves are smooth: `count_dissociation_tail_jumps` flags 0 pairs, against 7 PBE and 10 r2SCAN pairs before.
- **Models are not trained on symmetric dimer data.** Training sets contain almost no isolated dimers, and rattled bulk frames such as those of OMat24 have no symmetry to impose, so `ISYM` makes no difference to them. A model's dimer curve is an extrapolation either way, and the physical ground state is the fairer target for it.

Caveat: for about 30 pairs per functional the binding energy differs from the first version's `ISYM=2` curves by more than 50 meV (largest outside the broken Ho/Er/Tm curves: Dy, Pr, Sm, O, La, Pa, Pm, Zr), so diatomics metrics for those elements depend strongly on which of the two references they are scored against.

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

## Occupation Seeding

Within one spin candidate, an open shell still has many self-consistent states, one per pattern of occupied orbitals, and which one an SCF reaches depends on its start. For a free Ho atom the whole-electron 4f patterns span 111 meV in PBE and 379 meV in r2SCAN. Symmetry-broken SCF from a random start finds the PBE ground state of most f atoms to within 1 meV, but not reliably: Pm misses it by 42 meV in PBE, and Ho by 291 meV in r2SCAN.

Cold starts can therefore carry an occupation seed: each spin's shell occupation is written into the start's one-center occupation matrices as whole electrons in chosen real orbitals, held through the non-self-consistent delay steps, and then released. This is the occupation-matrix control of the DFT+U literature without the +U, so it biases only the start; every published energy is a free SCF.

- **f shells (14 open-f elements: Ce-Tm, Pa, U).** For each element and functional, every whole-electron f pattern of the free atom's ground configuration and its f^(n±1) neighbors was converged (666 free-atom calculations per functional). The lowest one or two distinct patterns seed chains on both atoms (the AFM partner with the spins swapped), swept in both directions. r2SCAN chains use the PBE patterns where those lie within 2 meV of r2SCAN's own lowest pattern, and r2SCAN's lowest pattern for Gd (f⁷ where PBE prefers f⁸), Tb, Pa, Dy and U.
- **p shells (18 open-p elements: B, Al, Ga, In, Tl; C, Si, Ge, Sn, Pb; O, S, Se, Te, Po; F, Cl, Br, I, At).** The partly filled spin's distinct orbital (the electron of p¹ and p⁴, the hole of p² and p⁵) is seeded either along the bond or perpendicular to it, both sweep directions, for every spin candidate except the closed-shell one. In Br₂ (PBE) the perpendicular hole lies up to 47 meV below every unseeded chain at 5 Å, and for the p¹ elements Al, Ga, In and Tl the seeded orientation lowers the dissociation tail by 16-46 meV. Half-filled p³ shells are spherical and need no seed.

## Sweeps

Every point warm-starts from its neighbor's wavefunction; the constant box keeps FFT grids compatible. Warm starts follow one electronic state, which is what makes curves smooth, but also means a chain keeps the state it started in. The reference therefore combines several kinds of chains per spin candidate:

- **Inward and outward sweeps.** One chain starts cold at 6 Å and sweeps inward, another starts cold at the bonded anchor (nearest `2 * r_cov`) and sweeps outward to 6 Å and inward to the wall. For the open-f elements the seeded chains run in both directions and the unseeded ones in one.
- **Inner sweeps.** A state that is lowest only at short range can escape both: the inward chain arrives in a long-range state and keeps it until that state gives way, and the outward chain starts beyond the short-range region. Where a merged curve showed an interior bump above 20 meV, a tail kink or a lower `ISYM=2` value, a chain per candidate starts cold four points inside the minimum and sweeps inward and then outward to 6 Å, and two more start two points to either side of each bump, in the state that wins there. In r2SCAN Ba₂ this follows the moment-free state from 3.2 Å outward and lowers the point at 3.4 Å, where the inward-sweeping chains had switched states, by 0.38 eV.
- **Rescue chains.** Where the first version's `ISYM=2` curve lies lower outside the repulsive wall, a chain starts cold at the distance of the largest gap, in that curve's spin candidate and site moments, and sweeps both ways.
- **Symmetric chains.** `ISYM=2` chains with otherwise identical settings, datasets and grid are merged too. They supply points only where they are lowest, which is mainly where breaking symmetry gains nothing (closed shells, where both agree to 0.1 meV).

Non-converged points never enter the merge.

Not every bump is an SCF artifact. PBE Al₂ has two near-degenerate triplet states with minima near 2.46 Å and 2.75 Å; where they cross, the merged curve has a 21 meV cusp at 2.65 Å. Seeded cold starts at 2.55 Å and 2.65 Å, in the orientation that finds the long-bond state further out, converge to the short-bond state there. The other bumps above 20 meV (PBE Pm₂, Hf₂; r2SCAN Zr₂, Eu₂, Tc₂, Th₂) likewise sit where the lowest state switches between spin candidates or orbital states; for Hf₂, Tc₂ and Zr₂, chains started cold on either side of the bump reproduce it.

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

The artifact contains 92 PBE and 92 r2SCAN homonuclear entries. Every distance has a converged energy except the `0.8 * r_cov` point of Eu and Gd, where the SCF diverges in both functionals (the first version has no energy there either). Merged curves can be shorter than 50 points because postprocessing drops 7 isolated points as likely SCF artifacts.

The quality checks are diagnostic and do not alter endpoints. `count_dissociation_tail_jumps` inspects adjacent energy steps among up to the final three merged points, flags steps of at least 0.1 eV, and ignores the short-range repulsive wall. `reference-quality.json` records the number of flagged steps as `tail_jumps`; the `tail_jump_pairs` summary counts each affected element-functional pair once. The artifact flags no pairs (the first version's `ISYM=2` curves flag 7 PBE and 10 r2SCAN pairs). Known residual gaps where the first version's `ISYM=2` curve found a lower state outside the repulsive wall: Sm/PBE by 163 meV at 4.7 Å, Fe/r2SCAN by 127 meV at 1.56 Å, and 6-69 meV at single distances in 8 other PBE and 9 other r2SCAN pairs (largest: Y/r2SCAN 69, Gd/r2SCAN 48, Ni/r2SCAN 40, La/PBE 39, Br/PBE 38 meV).

## Postprocessing

The postprocess is deliberately narrow. It applies local heuristics to points that are strongly suggestive of SCF artifacts while trying not to remove physical features such as barriers or state crossings. These classifications are not proofs: a sufficiently narrow physical feature could resemble a one-grid-point artifact. Anything that does not meet the explicit thresholds is published as computed, and the per-point projected moments and selected candidates expose residual branch sensitivity.

It handles three cases:

- candidate points that fall tens to thousands of eV below their neighboring branch, seen mostly for f-electron heavy lanthanides, are dropped as numerically implausible discontinuities before merging because they would otherwise always win the per-distance minimum;
- isolated severe (&ge;3 eV) drops in the merged curve are replaced by the neighboring candidate when the lower state appears at only one sampled separation; this treats a one-point notch as a likely branch failure
- isolated upward bumps (&ge;0.1 eV single-grid-point local maxima) are dropped as likely SCF artifacts, not because physical PECs are forbidden from having barriers or multiple wells.

The postprocess deliberately does not replace short, near-degenerate (&le;0.2 eV) "spin-branch islands" with the surrounding branch: that would edit far more points than all other filters combined for negligible smoothness gain, and would replace the computed minimum with a higher-energy branch.

Tests for this logic live in `tests/metrics/diatomics/test_diatomics_reference.py` and `tests/scripts/test_build_diatomic_reference.py`.

## Acknowledgements

Thanks to [Andrew S. Rosen](https://cbe.princeton.edu/people/andrew-rosen) [[Google Scholar](https://scholar.google.com/citations?user=lHBjgLsAAAAJ&hl=en)] for early advice on computing dimer energy curves and for suggesting the per-atom magnetic moment view.
