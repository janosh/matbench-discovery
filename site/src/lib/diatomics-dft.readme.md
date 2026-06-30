# Homonuclear Diatomic DFT Reference Curves

`diatomics-dft.json.gz` contains VASP potential-energy curves E(r) for homonuclear dimers X2 (X = H to U, Z=1-92) at PBE and r2SCAN, used as the DFT overlay on the Diatomics page. Format after decompression: `{ functional: { "<El>-<El>": { distances: number[], energies: number[], forces?: number[][][] } } }`; energies are sigma-to-0 (E0) values in eV.

The choices below were needed for smooth, physical curves; relaxing them caused cliffs, spikes, multi-eV jumps, or wrong dissociation limits.

## DFT setup

VASP 6 (GPU) inputs came from ferrox's MP24 static set: `ferrox vasp write-inputs --set mp24-static --xc-functional {pbe,r2scan} --potcar-mode real`. Per-point overrides:

- ISMEAR=0, SIGMA=0.05; KSPACING=0.6 (Gamma). A 2x2x2 check changed energies by <4 meV mean, max ~20 meV for diffuse atoms at large r.
- ISPIN=2 with pinned NUPDOWN/MAGMOM; NELM=250 for hard stretched/repulsive SCF.
- LWAVE=.TRUE. for warm starts; LCHARG/LAECHG/LVTOT/LELF=.FALSE. and LORBIT deleted to avoid ~180 MB per-point grid output.
- ISTART=1, ICHARG=0 with a seed wavefunction, else 0 and 2.

Box: 15 Å cubic vacuum. 13 Å left ~0.15 eV periodic-image error on diffuse atoms; 15 Å lowers this to ~0.03 eV; 18 Å was converged but ~2.7x costlier.

Distance grid: 50 geomspaced separations from 0.8\*r_cov to 6 Å. Going deeper caused severe PAW augmentation-sphere overlap and hundreds-of-eV repulsive energies, degrading the DFT reference.

## Warm-start: reuse the wavefunction from x+-dx for x

Cold SCF at large separation is near-degenerate open-shell atoms and barely converges. Instead, cold-start the bonded anchor (nearest 2\*r_cov), then warm-start outward and inward. Each point reuses `wavecar`/`vaspwave.h5` from its neighbor one step toward the anchor (`x+dx` inward, `x-dx` outward); the constant box keeps FFT grids compatible and reaches the hardest large-r SCF last.

## Spin: NUPDOWN pinning plus 2-candidate adiabatic

Cliffs and spikes mainly came from SCF hops between spin multiplicities. Each sweep pins total spin (`NUPDOWN=2S`) with ferromagnetic MAGMOM summing to NUPDOWN. Two fixed multiplicities are run per dimer and merged by per-point minimum:

- Molecular ground-state 2S (the experimental multiplicity) gives the correct well. Main group from MO theory; 3d metals from Gutsev and Bauschlicher, J. Phys. Chem. A 2003 (Cr2 singlet, Fe2 septet, Mn2 11-plet); known 4d/5d hardcoded (Mo2, W2 singlet).
- Ferromagnetic atomic 2S = 2\*(free-atom Hund moment) gives the correct open-shell dissociation; a single restricted multiplicity can dissociate to the wrong limit (restricted-singlet N2 is ~6 eV too high).
- Uncertain 4d/5d/4f/5f use a coarse {0, ~Hund, 2\*Hund} scan; Mn2 uses DFT high-spin (11-plet) because its antiferromagnetic ground state is not single-determinant accessible.

Antiferromagnetic MAGMOM init did not fix singlet dissociation because warm steps read the wavefunction and ignore MAGMOM; the FM candidate is needed.

## Reconciliation: escape warm-start branch traps

A post-sweep pass treats interior local maxima as warm-start branch traps (e.g. Cr2 orbital-occupation branches). It re-runs each bump seeded from its lower-energy neighbor, iterating until multi-point bumps collapse. To avoid f-electron thrashing: only bumps >= 0.3 eV are retried, unimproved points are not retried, with caps of 15 re-runs and 4 passes.

## Pipeline

1. pec*sweep.py: one (element, functional, NUPDOWN) sweep (anchor cold-start, bidirectional warm-start, reconciliation). Writes `<El>*<xc>\_n<NUPDOWN>/curve.json`.
2. gen_tasks.py: emits per-dimer candidate multiplicities and adiabatic_cands.json.
3. merge*adiabatic.py: per-point minimum over a dimer's candidates into the adiabatic ground-state `<El>*<xc>/curve.json`.
4. build_dft_references.py: aggregates merged curves into diatomics-dft.json.gz (sigma-to-0 energy per point, non-finite points dropped).
