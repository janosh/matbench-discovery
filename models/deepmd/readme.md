# DeePMD

DPA3 is an advanced interatomic potential leveraging the message-passing architecture, implemented within the DeePMD-kit framework, available on [GitHub](https://github.com/deepmodeling/deepmd-kit/tree/dpa3-alpha).

## DeePMD-DPA3

### Model

- DPA-3.1-3M-FT：model pretrained on the `OpenLAM` datasest (using DPA-3.1-3M from [DPA3-paper](https://arxiv.org/abs/2506.01686)), and finetuned with `MPtraj` & `Sub-Alex` datasets. Superseded by DPA-4.1-Max-OAM.

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

DPA4 is an SE(3)-equivariant interatomic-potential architecture built on an EMFA
(Edge-conditioned, Multi-Focus, Attention) SO(2)-equivariant convolution. It
combines a low-rank edge–node SO(2)-equivariant product, a multi-focus design
for message nonlinearity, and envelope-gated attention for message aggregation.
A Lebedev-grid projection preserves SO(3)-equivariance in the nonlinearity to
machine precision, and a compiler-friendly conservative energy-gradient training
path gives up to ~3× wall-clock speedup under `torch.compile`. See the
[DPA4 paper](https://arxiv.org/abs/2606.02419).

### DPA-4.0.1-Pro-MPtrj

DPA-4.0.1-Pro-MPtrj is the DPA4-Pro model trained only on MPtrj for this
submission; it supersedes DPA-4.0-Pro-MPtrj. Following the DPA4-Pro column of
Table S-17 in the paper, the model parameters are: feature dim 64, 6 interaction
layers, 2 focuses, 4 SO(2) layers, a rank-2 per-channel edge–node product, a
Bessel radial basis with 16 bases, lmax=5, mmax=1, Lebedev quadrature, SiLU+GLU
activations, a 6.0 Å cutoff, and up to 384 selected neighbors (~22.8M
parameters). It was trained for 2×10⁶ steps with the HybridMuon optimizer, a WSD
learning-rate schedule, MAE energy/force/virial loss weights 20/20/5, bf16 AMP,
TF32 matmul, and `torch.compile` enabled.

### DPA-4.1-OAM

DPA-4.1-Max-OAM, DPA-4.1-Air-OAM and DPA-4.1C-Air-OAM are trained on MPtrj,
OMat24 and sAlex (OSOD), with checkpoints released under CC-BY-4.0.
DPA-4.1-Max-OAM supersedes DPA-3.1-3M-FT. Max and Air
use the DPA4 architecture; C-Air uses the compact DPA4C architecture. Each model was trained
in two stages, pre-training and then fine-tuning, with the HybridMuon optimizer.
DPA-4.1-Max-OAM-PFT is the Max checkpoint after phonon fine-tuning. Step counts,
hardware and training costs are recorded in the model cards.

| Model | Parameters | Checkpoint and predictions |
| --- | ---: | --- |
| [DPA-4.1-Max-OAM](dpa-4.1-max-oam.yml) | 99,009,258 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Max-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755060) |
| [DPA-4.1-Air-OAM](dpa-4.1-air-oam.yml) | 5,147,909 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Air-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755042) |
| [DPA-4.1C-Air-OAM](dpa-4.1c-air-oam.yml) | 433,673 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1C-Air-OAM_model_checkpoint_and_Matbench_Discovery_predictions/34755132) |
| [DPA-4.1-Max-OAM-PFT](dpa-4.1-max-oam-pft.yml) | 99,005,146 | [Figshare](https://figshare.com/articles/dataset/DPA-4_1-Max-OAM-PFT_model_checkpoint_and_Matbench_Discovery_predictions/34846080) |

#### Inference

The released files are training checkpoints. The shared calculator exports each
checkpoint to a `.pt2` archive on first use (`dp --pt freeze` for DPA4,
`dp --pt-expt freeze` for DPA4C) and evaluates it with DeePMD-kit's `DeepPot`
and `DP.calculate`. Exports are cached by checkpoint content, backend, GPU model
and compute capability, and kernel policy. The first export of the Max model can
take a long time because of kernel tuning. Set `DP_TRITON_INFER` and
`DP_CUDA_INFER` before the export to choose the kernel policy, which is baked
into the archive. DeePMD-kit selects its device at import time: use
`CUDA_VISIBLE_DEVICES` to select one GPU, or set `DEVICE=cpu` for CPU inference.

DPA4 and DPA4C, including DPA-4.1-Max-OAM-PFT, use the uv project
`models/deepmd/dpa41-env` with PyTorch 2.13.0+cu132 and Python 3.14. DeePMD-kit
is installed from the `master` branch, so the environment follows the current
source rather than a fixed commit. The project builds DeePMD-kit from source
and needs Linux x86-64, a CUDA toolkit with `nvcc` and a C++ compiler. Run the shared runners as described in the contributing guide, e.g.

```sh
uv run models/run_discovery.py --print-cmd --model dpa-4.1-max-oam --dry-run
```
