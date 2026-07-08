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
    }
  }
}
```

`functional` is `PBE` or `r2SCAN`; `formula` is e.g. `O-O`. `distances`, `energies`, and `forces` share the same first axis over curve points. Energies are final sigma-to-0 VASP energies in eV; forces are per-point, per-atom Cartesian 3-vectors in eV/Å with axis order `[point][atom][x|y|z]`.

## DFT Setup

Each candidate curve was produced with VASP 6 using Materials Project (MP24)-style static input settings and real-space PAW potentials for PBE and r2SCAN.

Important per-point settings:

- `ISMEAR=0`, `SIGMA=0.05`.
- `KSPACING=0.6` (Gamma-only in this large cell). A 2x2x2 $k$-grid spot check on Na/PBE and Al/PBE changed final energies by 2.3 meV mean absolute over 40 matched points; Na's diffuse 6 Å tail (the largest separation) reached 23 meV, while Al stayed below 6.2 meV.
- `ISPIN=2` with fixed `NUPDOWN` for every spin candidate and explicit `MAGMOM` initialization where needed.
- `NELM=250` for difficult stretched or repulsive geometries.
- `LWAVE=.TRUE.` so neighboring geometries can reuse the wavefunction.
- `LCHARG`, `LAECHG`, `LVTOT`, and `LELF` disabled, and `LORBIT` removed, to avoid roughly 180 MB of unused grid output per point.
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

The artifact contains 92 PBE and 92 r2SCAN homonuclear entries. Known caveats: `Er/r2SCAN` has two short spin candidates (`NUPDOWN=2` with 43 finite points and AFM with 44 finite points, below the 45-point completeness threshold); `Ho/r2SCAN` and similar heavy lanthanides and actinides remain branch-trapped and jumpy despite all postprocessing and checks. These seem irrecoverable, at least with current VASP 6.4 pseudo-potentials.

## Postprocessing

The postprocess is deliberately narrow. It is meant to fix isolated SCF artifacts without smoothing away real physics like spin crossovers.

It currently handles three cases:

- isolated severe spin-branch drops, replacing the dropped point with the smoother neighboring branch;
- short, low-gain spin-branch islands, only when the surrounding branch is locally smoother than the adiabatic minimum point;
- isolated upward energy bumps.

Tests for this logic live in `tests/metrics/diatomics/test_diatomics_reference.py` and `tests/metrics/diatomics/test_build_diatomic_reference.py`.
