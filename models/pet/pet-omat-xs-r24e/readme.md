# PET-OMat-XS R24E

R24E shares the conservative PET-OMat-XS v1.0.0 weights across 24 orientations. It averages the energy over the proper tetrahedral group and its spatial inverses. The atomic forces and bulk stress are derivatives of this energy. The fitted energy correction adds one scalar per atom and one reference for its element. It has zero derivative with respect to coordinates and strain. The parent checkpoint is pinned in the model YAML. The elemental references and their training-data hashes are in [frozen-model.json](frozen-model.json).

UPET 0.3.1 already offers rotational averaging with `rotational_average_order` and a configurable batch size on a Lebedev–Laikov grid. R24E specifies a fixed 24-member group and adds composition references fitted on MPtrj. The comparison below uses the same checkpoint in its unaveraged conservative mode as the parent. Performance relative to the native UPET grid has not been measured. [UPET calculator](https://github.com/lab-cosmo/upet)

## Reproduce the fit

The 871 fit structures represent 438 reduced compositions. Each composition contributes the same total weight. The regression fits DFT energy per atom minus raw orientation-averaged PET energy per atom, using a Huber threshold of 0.1 eV/atom, six IRLS steps and a ridge penalty of 0.0005 on the element coefficients. The intercept is unpenalized. The package includes the atomic numbers, cell matrices, fractional coordinates and both total energies (eV) in the exact fitting order. The 48 independent confirmation compositions are excluded. Atom order, lattice vectors and fractional coordinates come from the MPtrj records without structural relaxation by either model.

From the repository root, fit from stored energies or recompute the orientation-averaged parent energies first.

```sh
uv run models/pet/train_pet_r24e.py
uv run models/pet/train_pet_r24e.py --recompute-parent --device cpu
```

Structure and DFT energy records are adapted from Bowen Deng's Materials Project Trajectory Dataset, version 2, Figshare file 41619375, published under the MIT license reproduced below. The per-frame identifiers correspond to the keys in the original archive. Raw PET energies use the parent checkpoint identified by SHA-256 in the frozen model. [MPtrj dataset and license](https://doi.org/10.6084/m9.figshare.23713842.v2) · [Parent model](https://huggingface.co/lab-cosmo/upet)

## Run the shared benchmarks

Use float32 to reproduce the reported model. These commands print the isolated commands for each shared runner. The printed command can be run locally or on a GPU worker with the same environment. Add `--device cuda` for an explicit GPU request. Discovery and phonons support resumable shards and strict merging through the [contribution guide](../../../contributing.md).

```sh
uv run models/run_discovery.py --model pet-omat-xs-r24e --dtype float32 --dry-run --print-cmd
uv run models/run_kappa.py --model pet-omat-xs-r24e --dtype float32 --dry-run --print-cmd
uv run models/run_md.py --model pet-omat-xs-r24e --dtype float32 --dry-run --print-cmd
uv run models/run_diatomics.py --model pet-omat-xs-r24e --dtype float32 --dry-run --print-cmd
```

## Development and pilot measurements

The 48 held-out MPtrj compositions each contributed two structures. Relative to the parent, energy MAE fell from 0.08694 to 0.04178 eV/atom and force MAE from 0.05750 to 0.05099 eV/Å. The force MAE ratio was 0.887 (paired-composition bootstrap 95% interval 0.838–0.920). The ratio for stress change was 0.993 (0.979–1.003).

On a fixed 50-structure WBM pilot, both models converged on all structures. Of 39 unique prototypes, five were DFT-stable. F1 was 0.286 for PET and 0.833 for R24E. The component ablation below uses the same 39 prototypes. The fitted element references account for most of the F1 and MAE gain on this pilot.

| Variant | F1 | Hull-distance MAE (eV/atom) | RMSE (eV/atom) |
| --- | ---: | ---: | ---: |
| Parent PET | 0.2857 | 0.11374 | 0.18698 |
| Parent plus elemental references | 0.8333 | 0.06977 | 0.13950 |
| Raw 24-orientation average | 0.2857 | 0.11841 | 0.19215 |
| R24E | 0.8333 | 0.06792 | 0.14159 |

On BaTe, BAs and InAs, mean κ-SRME decreased from 1.341 to 0.914. Both models substantially underestimated BAs conductivity. InAs had imaginary modes with the parent and a mode-resolved κ-SRME of 0.810 with R24E.

Two 0.25 ps NVT trajectories gave mixed observable changes. For anthracene the radial-distribution error increased from 10.236 to 11.102 and the angular-distribution error decreased from 0.235 to 0.231. For MAPbBr₃ the mean-pressure error decreased from 0.378 to 0.369 GPa, while the pressure-distribution Wasserstein distance increased from 0.436 to 0.443 GPa. R24E took 12.7 and 15.4 times the parent CPU runtime for these two systems. The full task calls for 17 trajectories of 20 ps. The official model metadata is pending.

## Data license

MPtrj records are attributed to Bowen Deng, 2023, Materials Project Trajectory Dataset v2. The following MIT grant applies to the redistributed structure and DFT-energy records.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
software.
