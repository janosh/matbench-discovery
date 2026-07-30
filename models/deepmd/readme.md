# DeePMD

DPA3 is an advanced interatomic potential leveraging the message-passing architecture, implemented within the DeePMD-kit framework, available on [GitHub](https://github.com/deepmodeling/deepmd-kit/tree/dpa3-alpha).

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
pip install git+https://github.com/deepmodeling/deepmd-kit@devel
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
