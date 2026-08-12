# EMA-GNN

EMA-GNN is a six-model structural GNN ensemble that predicts relaxed formation energy per atom directly from unrelaxed structures (`IS2E`, target `E`), without relaxation at inference time.

## Architecture and training

- 2,209,793 parameters per estimator, three message-passing blocks, 256 hidden features, and a 4.0 Å graph radius
- Sequential edge, node, and global updates with adjacency-normalized edge-to-node messages
- Six independently seeded estimators trained for 500 epochs on Materials Project 2022 relaxed structures and evaluated with EMA weights

The architecture follows the structural GNN in Merchant et al., ["Scaling deep learning for materials discovery"](https://www.nature.com/articles/s41586-023-06735-9). The existing GNoME leaderboard entry is the paper's separate NequIP-type interatomic potential. EMA-GNN uses six estimators instead of ten and 500 epochs instead of 1000. Exact hyperparameters, validation results, and ablations are recorded in [`ema-gnn.yml`](ema-gnn.yml).

Per-seed held-out MAEs were 23.32, 24.02, 24.09, 24.33, 23.78, and 24.41 meV/atom.

## Inference

Each estimator takes the minimum prediction across 20 isotropic volume augmentations linearly spaced from 0.8–1.2x volume; the ensemble combines estimators by median.

## Artifacts and reproduction

- [Model metadata, metrics, and predictions](ema-gnn.yml)
- [Checkpoints](https://doi.org/10.6084/m9.figshare.33111509)
- [Training and evaluation code](https://github.com/submerged-in-matrix/gnome-repro-structural)

Checkpoint loading requires `strict=False`: `avg_adjacency` is registered as a buffer but supplied through the constructor from the checkpoint's `stats` dictionary, so it is absent from `model_state`.

EMA-GNN is an archived direct predictor rather than an ASE calculator. `models/run_discovery.py` derives formation energy from calculator total energy after Materials Project compatibility corrections; passing EMA-GNN's formation-energy output through that route would apply incompatible corrections.
