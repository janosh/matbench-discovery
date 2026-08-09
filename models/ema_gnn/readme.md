# EMA-GNN

Structural graph neural network for direct formation-energy prediction from
unrelaxed crystal structures.

## Model

EMA-GNN maps an unrelaxed input structure to the formation energy per atom of the
corresponding relaxed crystal, with no relaxation step at inference time. Each
message-passing block updates edge, node and global features in sequence, with
edge-to-node messages normalized by the dataset-average adjacency so that message
magnitude does not scale with coordination number.

| Property | Value |
| --- | --- |
| Task | `IS2E` |
| Targets | `E` |
| Training set | Materials Project 2022, relaxed structures |
| Parameters | 2,209,793 |
| Ensemble | 6 independently seeded models |
| Cutoff radius | 4.0 Å |
| Hidden dimension | 256 |
| Message-passing blocks | 3 |

The architecture follows the structural GNN described in Merchant et al.,
["Scaling deep learning for materials discovery"](https://www.nature.com/articles/s41586-023-06735-9),
Nature 624 (2023). That paper describes two models; the entry on this leaderboard
under the name GNoME is the other one — a NequIP-type interatomic potential
submitted as `IS2RE-SR` with targets `EF_G`. The structural GNN benchmarked here
predicts energy directly from an unrelaxed input and had not previously been
submitted.

The published architecture was the starting point, not an assumption. An anchored
coordinate search over hidden-layer width, activation, network depth and learning
rate was run on a stratified subset held out from training, with WBM withheld
throughout. Every axis converged to the configuration reported in the paper.

## Training

Six models trained serially on a single NVIDIA RTX 4070 Ti, differing only in
random seed and sharing the same training split. Held-out MAE per seed, meV/atom:
23.32, 24.02, 24.09, 24.33, 23.78, 24.41 — a spread of 1.1 meV with no diverged
run.

```txt
optimizer          Adam
learning rate      5.5e-4, LinearLR to 0.1x final
epochs             500
batch size         128 x 2 accumulation steps (effective 256)
gradient clip      1.0
loss               L1 on standardized targets
EMA decay          0.999
early stopping     disabled
weights evaluated  EMA
```

Two deviations from the reference description, both compute-driven: six ensemble
members rather than ten, and 500 training epochs rather than 1000. Exponential
moving average of the weights is an addition to the published recipe; a single-seed
comparison recorded a best EMA validation MAE of 23.77 meV/atom at epoch 386 against
a final-epoch raw MAE of 24.31, and EMA continued improving past the point where the
raw weights plateaued, which is why early stopping was disabled.

## Inference

Twenty-point isotropic volume test-time augmentation per structure (0.8–1.2x
volume, linear in volume), aggregated by minimum per model, then by median across
the six models. This is the paper's protocol. Median rather than mean across
ensemble members prevents a single out-of-distribution failure from displacing the
ensemble prediction.

## Predictions

Predictions on the WBM test set are hosted on Figshare and linked from
[`ema-gnn.yml`](ema-gnn.yml). Checkpoints for all six seeds are archived in the
same item: <https://doi.org/10.6084/m9.figshare.33111509>

Loading a checkpoint requires `strict=False`, because `avg_adjacency` is registered
as a buffer but supplied through the constructor from the checkpoint's `stats` dict,
and is therefore absent from `model_state`.

## Notes for reviewers

EMA-GNN is an archived direct-prediction submission and does not expose an ASE
calculator for `models/run_discovery.py`. That runner obtains a total energy from a
calculator, applies `MaterialsProject2020Compatibility`, then subtracts Materials
Project elemental references to derive `e_form_per_atom`. EMA-GNN emits
`e_form_per_atom` directly, so routing it through the runner would require
inverting the reference subtraction and would apply MP2020 corrections on top of a
quantity that already implicitly contains them. Predictions are submitted as a
hosted artifact, following the precedent of other `IS2E` entries.

Training and evaluation code, including the scoring script used to produce the
metrics in `ema-gnn.yml`, is at
<https://github.com/submerged-in-matrix/gnome-repro-structural>.
