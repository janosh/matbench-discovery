# Möbius OAM

Möbius OAM is an equivariant atomistic model released by Zhongguancun Academy
and Zhongguancun Institute of Artificial Intelligence. The Matbench Discovery
entry uses the public `eq_0.1` inference checkpoint and keeps the benchmark key
`mobius-oam`.

The model is pretrained on OMat24 and fine-tuned on sAlex and MPTrj. It predicts
energies, conservative forces, and virials for its supported element set. The
public checkpoint and benchmark artifacts are available from the
[Figshare record](https://doi.org/10.6084/m9.figshare.35010493).

## Inference

Install the public [RUNE](https://github.com/Jonsnow-willow/rune) release and
use the ASE interface:

```python
from ase.io import read
from rune.mobius.plugin import build_calculator

atoms = read("structure.extxyz")
atoms.calc = build_calculator("model.pt", device="cpu")
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()
```

The pinned RUNE release `7f2a36c643bb79149b4714c04263002b7aea5f92` preserves the
dtype stored in the Möbius checkpoint. The public `model.pt` contains `float32`
weights, so this calculator runs in `float32`. This release does not expose a
Möbius `dtype` override; precision is checkpoint-defined. The historical
`float32` cost label records the actual inference precision and is not a
separate loader switch.

The exact checkpoint is linked by `models/mobius/mobius-oam.yml`. The RUNE code
is released under [LGPL-3.0](https://github.com/Jonsnow-willow/rune/blob/main/LICENSE),
and the inference checkpoint is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode).

## Benchmark protocol

The discovery and geometry-optimization artifacts were generated with Matbench
Discovery commit `7a7a941` using float32 inference, FIRE with
`FrechetCellFilter`, `max_force=0.02` eV/Å, and at most 500 optimization steps.
The replacement κ run uses current main commit
`6d587c184cdb4a4865dc0b15b04a5e6d77ed464f`, a 0.03 Å finite-displacement
distance, 300 K, `force_max=1e-4`, `max_steps=300`, and `save_forces=true`.
It includes the required harmonic phonon sidecar with schema v3 and the HDF5
force sidecar. The four κ artifacts are version 2 of the Figshare record.

## Results

| Task | Result |
| --- | ---: |
| Discovery F1 (full test set) | 0.866356 |
| Discovery F1 (unique prototypes) | 0.885614 |
| Geometry RMSD | 0.062821 |
| κ<sub>SRME</sub> | 0.296886 |
| κ<sub>SRE</sub> | 0.191703 |
| κ<sub>SRD</sub> | -0.066578 |
| κ failure rate | 0.009709 |
| Imaginary-mode rate | 0.009709 |

Contact: Jiahui Liu, Zhongguancun Academy and Zhongguancun Institute of Artificial Intelligence,
[liujiahui@zgci.ac.cn](mailto:liujiahui@zgci.ac.cn).
