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

The exact checkpoint is linked by `models/mobius/mobius-oam.yml`. The RUNE code
is released under [LGPL-3.0](https://github.com/Jonsnow-willow/rune/blob/main/LICENSE),
and the inference checkpoint is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode).

## Historical benchmark protocol

The published artifacts were generated with Matbench Discovery commit `7a7a941`
using float32 inference, FIRE with `FrechetCellFilter`, `max_force=0.02` eV/Å,
and at most 500 optimization steps. The κ calculation used a 0.01 Å
finite-displacement distance, 300 K, `force_max=1e-4`, and at most 300 steps.
The harmonic sidecar is not part of this submission because the recovered file
uses an older schema than the current Matbench loader.

## Results

| Task | Result |
|---|---:|
| Discovery F1 (full test set) | 0.866356 |
| Discovery F1 (unique prototypes) | 0.885614 |
| Geometry RMSD | 0.062821 |
| κ<sub>SRME</sub> | 0.299675 |

Contact: Jiahui Liu, Zhongguancun Academy and Zhongguancun Institute of Artificial Intelligence,
[liujiahui@zgci.ac.cn](mailto:liujiahui@zgci.ac.cn).
