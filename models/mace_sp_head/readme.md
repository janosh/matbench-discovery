# MACE-MP-0 SP-Head

MACE-MP-0 SP-Head predicts the relaxed formation energy of a WBM crystal from its **unrelaxed**
structure. It needs one frozen MACE-MP-0 (MPtrj-only) single point per structure and no geometry
optimisation. A LightGBM head, trained on MP 2022 data only, adds the missing relaxation energy.
It works from the single-point forces, stress and invariant MACE node features. The training set
is relaxed MP structures plus perturbed copies, all labelled with the relaxed DFT energy.

Unique prototypes, 3-seed mean ensemble:

| model | F1 | DAF | MAE (eV/atom) | R² |
|---|---|---|---|---|
| this model | 0.654 | 4.344 | 0.074 | 0.392 |
| raw MACE-MP-0 single point, no head | 0.554 | 4.699 | 0.129 | — |

Per-seed F1 is 0.6509, 0.6544 and 0.6526.

Code, trained heads, research log and reproduction steps:
<https://github.com/Ephemeral6/mace-sp-head>

- `train_mace_sp_head.py` builds the MP 2022 training set, runs the MACE single points and trains
  the LightGBM heads.
- `test_mace_sp_head_discovery.py` produces the WBM predictions.

Built with AI assistance: an Anthropic Claude agent designed, ran and documented the method.
