# BAM-MP-DeNS

BAM-MP-DeNS is an equivariant message-passing interatomic potential built with the
RACE architecture in [BAM-torch](https://github.com/myung-group/BAM-torch) (branch
`BAM_DeNS`). It was trained on the MPtrj dataset with denoising non-equilibrium
structures (DeNS) as an auxiliary task and predicts energy, forces, and stress (forces
and stress from autograd of the energy).

## Architecture

- RACE equivariant GNN, 5 message-passing layers, `max_ell=3`
- Hidden features `128x0e + 128x1o + 128x2e + 128x3o`, 12 radial basis functions,
  radial polynomial envelope `p=6`
- Cutoff 7.0 Å, `avg_num_neighbors=92`, 89 species, periodic boundary conditions
- 6,468,994 trainable parameters — a single model, not an ensemble or multi-head

## Training data and splits

- Trained **only** on MPtrj (Materials Project relaxation trajectories).
- MPtrj was split into separate training and validation sets (1,580,395 frames in
  total); no other data source was used.

## Training configuration

- Training ran for 118 epochs in restarted segments (`restart: true`) on 4 NVIDIA
  A100 40GB GPUs, with DeNS enabled throughout (probability 0.4, noise sigma uniform in
  [0.01, 0.1]). The submitted checkpoint is the best-validation checkpoint at epoch 118.
- Loss: Huber (energy `delta=0.05`, force `delta=0.25`, stress `delta=0.05`) with
  weights energy 1, force 1, stress 1
- Optimizer Adam (learning rate 1e-3) with ReduceLROnPlateau, EMA decay 0.999

## Checkpoint

`BAM-MP-DeNS-model.pkl` (156 MB, md5 `06f303110248a263f9ea2f4d2775af93`) is a full
training checkpoint (model weights, Adam optimizer state, scheduler, and EMA state).
Inference uses the single 6.5M-parameter model with its EMA weights, so
`model_params` in the YAML reports 6,468,994.

Loading this checkpoint requires `bam-torch` at BAM_DeNS commit `24ed0f0`, which passes
the DeNS-era architecture options from the checkpoint config to the model.

The benchmark predictions were computed with Python 3.13. The YAML sets
`python_version: '3.14'` because the shared runners require Python >= 3.14; all four
runners were smoke-tested with it.

## No-leakage statement

BAM-MP-DeNS was trained and selected **solely on the MPtrj train/validation splits**
(every training frame carries a Materials Project material ID). No WBM, convex-hull,
or Matbench Discovery test data was used at any stage of training, normalization, or
hyperparameter selection.

## Relaxations

Geometry optimizations used ASE `FIRE` with `FrechetCellFilter`, `max_force=0.05`
eV/Å or 500 steps.
