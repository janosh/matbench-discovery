# BAM-MP-core

BAM-MP-core is an equivariant message-passing interatomic potential built with the
RACE architecture in [BAM-torch](https://github.com/myung-group/BAM-torch). It was
trained from scratch on the MPtrj dataset for Matbench Discovery and predicts energy,
forces, and stress (`free_energy` is set equal to `energy` for ASE compatibility).

## Architecture

- RACE equivariant GNN, 5 message-passing layers, `max_ell=3`
- Hidden features `128x0e + 64x1o + 32x2e + 32x3o`, 8 radial basis functions
- Cutoff 6.0 Å, `avg_num_neighbors=60`, 89 species, periodic boundary conditions
- 17,023,365 trainable parameters — a single model, not an ensemble or multi-head

## Training data and splits

- Trained **only** on MPtrj (Materials Project relaxation trajectories).
- MPtrj was split into separate training and validation sets; no other data source was used.
- The per-element reference energies (`enr_avg_per_element`) and the energy
  scale/shift terms (`train_scale_shift`, `valid_scale_shift`) were fitted
  **exclusively on the MPtrj train/validation splits**.

## Training configuration

- 48 epochs, trained on NVIDIA H100 GPUs
- Loss: Huber (`delta=0.01`) with weights energy 10, force 1, stress 10

## Checkpoint

`BAM-MP-core-model.pkl` (409 MB, md5 `0d7182820d24ab51a452d475c9e5001d`) is a full
training checkpoint. Its size reflects bundled training state, not additional
parameters:

| component | size |
| --- | --- |
| model weights + buffers (fp32) | ~69 MB |
| Adam optimizer state | ~204 MB |
| EMA shadow + collected params | ~136 MB |

Inference uses the single 17M-parameter model with its EMA weights, so
`model_params` in the YAML reports 17,023,365.

## No-leakage statement

BAM-MP-core and all of its fitted normalization terms (per-element reference
energies and energy scale/shift) were derived **solely from the MPtrj
train/validation splits**. No WBM, convex-hull, or Matbench Discovery test data was
used at any stage of training, normalization, or hyperparameter selection.

## Relaxations

Geometry optimizations used ASE `FIRE` with `max_force=0.05` eV/Å or 500 steps.
