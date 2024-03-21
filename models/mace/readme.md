## MACE formation energy predictions on WBM test set

The original MACE submission used the 2M parameter checkpoint [`2023-08-14-mace-yuan-trained-mptrj-04.model`](https://figshare.com/ndownloader/files/42374049) trained by Yuan Chiang on the [MPtrj dataset](https://figshare.com/articles/dataset/23713842).
We initially tested the `2023-07-14-mace-universal-2-big-128-6.model` checkpoint trained on the much smaller [original M3GNet training set](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599) which we received directly from Ilyes Batatia. MPtrj-trained MACE performed better and was used for the Matbench Discovery submission.

In late October (received 2023-10-29), Philipp Benner trained a much larger 16M parameter MACE for over 100 epochs in MPtrj which achieved an (at the time SOTA) F1 score of 0.64 and DAF of 3.13.

### Convergence criteria

MACE relaxed each test set structure until the maximum force in the training set dropped below 0.05 eV/Å or 500 optimization steps were reached, whichever occurred first.

### Hyperparameters

#### Equivariance

- `hidden_irreps="64x0e + 64x1o + 64x2e"`
- `max_ell=3`
- `correlation=3`
- `num_interactions=2`
- `num_cutoff_basis=10`

#### Training

See the module doc string in `train_mace.py` for how to install MACE for multi-GPU training.
A single-GPU training script that works with the current [MACE PyPI release](https://pypi.org/project/mace-torch) (v0.3.4 as of 2024-03-21) could be provided if there's interest.

Our training used conditional loss weighting. We did _not_ use MACE's newest attention block feature which in our testing performed significantly worse than `RealAgnosticResidualInteractionBlock`.
