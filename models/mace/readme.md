## MACE formation energy predictions on WBM test set

This submission uses the [`2023-08-14-mace-yuan-trained-mptrj-04.model`](https://figshare.com/ndownloader/files/42374049) checkpoint trained by Yuan Chiang on the [MPtrj dataset](https://figshare.com/articles/dataset/23713842).
We initially tested the `2023-07-14-mace-universal-2-big-128-6.model` checkpoint trained on the much smaller [original M3GNet training set](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599) which we received directly from Ilyes Batatia. MPtrj-trained MACE performed better and was used for the Matbench Discovery v1 submission.

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

- `loss="uip"`
- `energy_weight=1`
- `forces_weight=1`
- `stress_weight=0.01`
- `r_max=6.0`
- `lr=0.005`
- `batch_size=10`

We used conditional loss weighting. We did _not_ use MACE's newest attention block feature which in our testing performed significantly worse than `RealAgnosticResidualInteractionBlock`.
