## MACE formation energy predictions on WBM test set

The original MACE submission used the 2M parameter checkpoint [`2023-08-14-mace-yuan-trained-mptrj-04.model`](https://figshare.com/ndownloader/files/42374049) trained by Yuan Chiang on the [MPtrj dataset](https://figshare.com/articles/dataset/23713842).

We also tested the `2023-07-14-mace-universal-2-big-128-6.model` checkpoint trained by Ilyes Batatia on the [MPF training set](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599). Although both derived from the Materials Project, MPF has both less materials and less frames than MPtrj. The MPtrj-trained MACE performed better and was used for the original Matbench Discovery submission.

On 2023-10-29 Philipp Benner shared a much larger 16M parameter MACE he trained for over 100 epochs in MPtrj which achieved an (at the time SOTA) F1 score of 0.64 and DAF of 3.13.

As of 2024-02-06 we have used the 4.7M parameter checkpoint [`2023-12-03-mace-128-L1.model`](https://figshare.com/ndownloader/files/42374052) trained by Yuan Chiang on the [MPtrj dataset](https://figshare.com/articles/dataset/23713842) for the leaderboard. This model was trained for the `MACE-MP-0` collaboration ([see manuscript](https://arxiv.org/abs/2401.00096)).

#### Training

The current "2023-12-03-mace-128-L1" model was trained using the MACE CLI with the training command recorded [here](https://github.com/ACEsuit/mace-mp/blob/main/mace_mp_0/2023-12-03-mace-128-L1.sh) in the `mace_mp_0` [repository](https://github.com/ACEsuit/mace-mp/tree/main).

For historical reasons an example training script is provided in `train_mace.py` showing how to train MACE with multiple gpus. Earlier iterations of this script were used to train the `2023-08-14-mace-yuan-trained-mptrj-04` and `2023-07-14-mace-universal-2-big-128-6` models. This script is not guaranteed to work with the current MACE PyPI release.

### Convergence criteria for UIP (Universal Interatomic Potential) relaxations

MACE relaxed each test set structure until the maximum force in the training set dropped below 0.05 eV/Å or 500 optimization steps were reached, whichever occurred first. The `FIRE` optimizer from `ase` was used with the `FrechetCellFilter` from `ase` to update the atomic positions.
