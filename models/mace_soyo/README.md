# MACE-SOYO

[MACE-SOYO](https://github.com/PrincessHakuryu/MACE-SOYO) is a multi-head interatomic potential built on the [MACE-MH-1](https://github.com/ACEsuit/mace) architecture. It adds attention-weighted edge weights, an equivariant RMSNorm, a heavier readout, the Muon optimizer and DPA4's warmup–stable–decay learning-rate schedule. Training uses `torch.compile` and inference an AOTInductor (AOTI) package, both CUDA-only (cuEquivariance, nvalchemi-toolkit). Installation, training, AOTI export, TorchSim and LAMMPS instructions live in the upstream repo.

## Checkpoints

Both v0.1 checkpoints are MIT-licensed [release assets](https://github.com/PrincessHakuryu/MACE-SOYO/releases/tag/v0.1):

- [`mace_soyo_mptraj.pth`](https://github.com/PrincessHakuryu/MACE-SOYO/releases/download/v0.1/mace_soyo_mptraj.pth) (MACE-SOYO-MPtrj): one `mptraj` head, trained on MPtrj.
- [`mace_soyo_oam.pth`](https://github.com/PrincessHakuryu/MACE-SOYO/releases/download/v0.1/mace_soyo_oam.pth) (MACE-SOYO-OAM): trained on OMat24, sAlex and MPtrj. Benchmarks use its `salexmp` head.

The shared runners (e.g. `uv run models/run_discovery.py --model mace_soyo_oam`) download the `.pth` checkpoint and export it once per GPU architecture to an AOTI `.pt2` package, which needs `nvcc`. The kappa and diatomics runners also accept an already exported `.pt2` via `--checkpoint`.

## Comparison with MACE-MH-1 on MPtrj

Both models were trained by the authors on the same MPtrj train/validation split for 100 epochs with E:F:S loss weights 1:1:1. Discovery metrics use unique prototypes; RMSD uses `symprec=1e-2`.

| Metric | MACE-MH-1 | MACE-SOYO-MPtrj |
| --- | ---: | ---: |
| CPS | 0.735 | 0.763 |
| F1 | 0.7549 | 0.7766 |
| DAF | 4.5574 | 4.7949 |
| Discovery MAE (eV/atom) | 0.04449 | 0.04136 |
| κ_SRME | 0.4197 | 0.3496 |
| RMSD | 0.0880 | 0.0832 |

On 4 A100s, MPtrj training took 5.0 GPU-days with compiled training (DPA4-Mini 5.8, DPA4-Neo 8.6, MACE-MH-1 with cuEquivariance 8.3). On a 576-atom NVT benchmark, an ASE step takes 11.2 ms, vs 22.7 ms for MACE-MH-1 with cuEquivariance and `torch.compile`.

Discovery and phonons ran on an RTX 4090, MD and diatomics on an H200 NVL.
