# MACE

MACE is a higher-order equivariant message-passing neural network for fast and accurate force fields.

Leaderboard checkpoints are `mace-mp-0` ([MACE-MP-0 paper](https://arxiv.org/abs/2401.00096),
[`2023-12-03-mace-128-L1.model`](https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0/2023-12-03-mace-128-L1_epoch-199.model))
and `mace-mpa-0`. The current MACE-MP-0 checkpoint was trained with
[this command](https://github.com/ACEsuit/mace-foundations/blob/f59cbe3e51/mace_mp_0/2023-12-03-mace-128-L1.sh).

The original submission used Yuan Chiang's 2M-parameter
[`2023-08-14-mace-yuan-trained-mptrj-04.model`](https://figshare.com/files/42374049),
trained on MPtrj. We also tested Ilyes Batatia's
`2023-07-14-mace-universal-2-big-128-6.model`, trained on
[MPF](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599), but the
MPtrj-trained model performed better. Philipp Benner later shared a 16M-parameter MACE
trained for over 100 MPtrj epochs, which achieved F1 0.64 and DAF 3.13.

`train_mace.py` is a historical multi-GPU example, not guaranteed to match current
MACE releases. Relaxations used ASE `FIRE` + `FrechetCellFilter` with
`max_force=0.05` eV/Å or 500 steps.
