# DeePMD

DPA3, DPA4 and DPA4C are interatomic potentials implemented in
[DeePMD-kit](https://github.com/deepmodeling/deepmd-kit).

## DeePMD-DPA3

### Model

- DPA-3.1-3M-FT：model pretrained on the `OpenLAM` datasest (using DPA-3.1-3M from [DPA3-paper](https://arxiv.org/abs/2506.01686)), and finetuned with `MPtraj` & `Sub-Alex` datasets.

  ```bash
  wget https://figshare.com/files/55141895
  ```

- DPA-3.1-MPtrj：model trained only on the `MPtraj` dataset.

  ```bash
  wget https://figshare.com/files/55141124
  ```

#### How to install

```bash
uv pip install git+https://github.com/deepmodeling/deepmd-kit@devel
```

#### How to use

```py
from ase import Atoms
from deepmd.calculator import DP

water = Atoms(
    "H2O",
    positions=[(0.7601, 1.9270, 1), (1.9575, 1, 1), (1.0, 1.0, 1.0)],
    cell=[100, 100, 100],
    calculator=DP(model="dpa-3.1-mptrj.pth"),
)

print(water.get_potential_energy())
print(water.get_forces())
```

### Results

Discovery metrics are recorded in each model YAML. Relaxed structures are available
from Figshare:

```sh
# DPA-3.1-3M-FT
wget https://figshare.com/files/55141109

# DPA-3.1-MPtrj
wget https://figshare.com/files/55141127
```

## DeePMD-DPA4

DPA4 is an SE(3)-equivariant interatomic-potential architecture that learns atomic
energies from chemical identities and three-dimensional atomic environments.
Its SO(2) convolutions couple edge and atom features across angular degrees,
combining a multi-focus representation with attention-weighted neighbor
aggregation. Wigner bilinear nonlinearities capture directional correlations
while preserving rotational symmetry. Forces and stresses are derived from the
predicted energy, giving a conservative potential for structure relaxation and
molecular dynamics. See the [DPA4 paper](https://arxiv.org/abs/2606.02419).

DPA4C is a compact, O(3)-equivariant architecture designed for efficient molecular
dynamics. A single local aggregation builds directional atomic features, which
are converted into invariant descriptors and mapped to atomic energies.
Separating radial functions from chemical-type dependence allows tabulated
compression for fast evaluation. Its strictly local representation avoids
exchanging learned atomic features between simulation domains, supporting
large-scale parallel simulations. Forces and stresses are derived from the same
energy. See the [DPA4C paper](https://arxiv.org/abs/2608.19041).

### DPA-4.0.1-Pro-MPtrj

DPA-4.0.1-Pro-MPtrj is a 22.8-million-parameter DPA4 model trained on MPtrj.
It supersedes DPA-4.0-Pro-MPtrj.

### DPA-4.1-OAM

The Max and Air models use DPA4; C-Air uses DPA4C. All three are trained on
MPtrj, OMat24 and sAlex.
DPA-4.1-Max-OAM-PFT adds phonon fine-tuning on Hessian labels from MDR-MP PBE ω_q,
with MPtrj co-training. All four checkpoints are released under CC-BY-4.0.
Model architectures use the `v20260911` presets specified in the model cards.

| Model | Parameters | Checkpoint and predictions |
| --- | ---: | --- |
| [DPA-4.1-Max-OAM](dpa-4.1-max-oam.yml) | 99,009,258 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Max-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755060) |
| [DPA-4.1-Air-OAM](dpa-4.1-air-oam.yml) | 5,147,909 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Air-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755042) |
| [DPA-4.1C-Air-OAM](dpa-4.1c-air-oam.yml) | 433,673 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1C-Air-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755132) |
| [DPA-4.1-Max-OAM-PFT](dpa-4.1-max-oam-pft.yml) | 99,005,146 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Max-OAM-PFT_model_checkpoint_and_Matbench_Discovery_predictions/34846080) |

Training costs in the model cards include all training stages, including
base-model training for PFT. The cards also report evaluation settings and limitations.
All evaluations used NVIDIA RTX PRO 6000 Blackwell Server Edition GPUs.

#### Inference

For ASE, load a checkpoint with `DP(model="<checkpoint>.pt", freeze=True)` to export and compile it for fast inference. For LAMMPS, see the [DPA4](https://github.com/deepmodeling/deepmd-kit/blob/master/doc/model/dpa4.md#export-and-running-in-lammps) and [DPA4C](https://github.com/deepmodeling/deepmd-kit/blob/master/doc/model/dpa4c.md#model-compression) docs. The benchmark runners use the environment in `models/deepmd/dpa41-env` (Linux x86-64, CUDA toolkit with `nvcc`), e.g. `uv run models/run_discovery.py --model dpa-4.1-max-oam`.
