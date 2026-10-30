# ALIGNN-FF

Atomistic Line Graph Neural Network force fields. Two entries live in this
family: the active `alignn-ff-mptrj` submission described first, and the
historical record of the aborted 2023 `alignn-ff` attempt kept below.

## `alignn-ff-mptrj` (active)

ALIGNN 2.0 universal force field
([arXiv:2609.19487](https://arxiv.org/abs/2609.19487)): a pure-PyTorch
reimplementation of ALIGNN that removes the DGL dependency. Messages pass
alternately over an atom graph and a line graph whose nodes are bonds and
whose edges carry bond angles. Energy is the only directly predicted quantity.
Forces and stress follow from one backward pass through the interatomic
displacement vectors, so the forces are conservative and the stress is the
analytic virial of the same energy.

- Checkpoint: `alignn_ff_mptrj.zip` in figshare item
  [10.6084/m9.figshare.34846446](https://doi.org/10.6084/m9.figshare.34846446)
  (CC BY 4.0), 4.06M parameters. Radius graph with a 5.0 Å cutoff, 3.5 Å
  three-body cutoff and at most 52 neighbors, 4 ALIGNN + 4 GCN layers at
  hidden width 256. Bond features are scaled by a smooth polynomial envelope
  that reaches zero at 4 Å.
- Per-element embedding: each atom gets its embedded CGCNN features plus a
  learned per-element vector. The CGCNN features alone are identical for Gd
  and Sm, Th and Np, La and Ce, and Ho and Er, so a model without the
  embedding cannot tell those elements apart. Loading the checkpoint needs an
  alignn release that includes the embedding, which `alignn-ff-mptrj.yml`
  pins.
- Training: from scratch on MPtrj with uncorrected energies, with energy per
  atom, forces and stress all in the L1 loss (weights 1, 1 and 0.1).
  1,422,355 frames for training and 79,019 each for validation and testing,
  30 epochs on 8 NVIDIA GB10 GPUs (effective batch 16, OneCycle with peak
  learning rate 3e-4). The submitted checkpoint has the lowest validation
  loss.
- Known limitations: MPtrj is PBE without dispersion, so layered van der
  Waals systems are poorly described.
- Precision: the calculator honors the shared runners' `--dtype` flag by
  setting torch's default dtype before the weights load, the same mechanism
  MACE uses. The weights were trained in float32.

The calculator is registered as `alignn_ff_mptrj` in
`matbench_discovery/calculators.py`. Discovery and geometry optimization use
the shared runner `models/run_discovery.py` with the isolated environment
declared in `alignn-ff-mptrj.yml`. Relaxations use ASE `FIRE` +
`FrechetCellFilter` with `max_force=0.05` eV/Å or 500 steps. Phonons use the
shared runner `models/run_kappa.py` with the `phonondb-v1` protocol. All
submitted predictions come from these runners at commit `c614eb0` and are
scored by the current evaluation code. 20 of the 103 materials have imaginary
modes after relaxation, so their conductivity is skipped and they count as
failures.

## `alignn-ff` (aborted, 2023, historical record)

The [ALIGNN FF model submission](https://github.com/janosh/matbench-discovery/pull/47) intended to get a complete set of formation energy predictions for the WBM test set post-ALIGNN-FF structure relaxation (i.e. the WBM IS2RE task).

This effort was aborted for the following reasons:

1. **Incompatibility issues**: ALIGNN-FF was pre-trained on the JARVIS data, which among other differences uses the OptB88vdW functional and is incompatible with the WBM test set generated using Materials Project workflows.
1. **Training difficulties**: ALIGNN-FF proved to be very resource-hungry. [12 GB of MPtrj training data](https://figshare.com/articles/dataset/23713842) turned into 600 GB of ALIGNN graph data. This forces small batch size even on nodes with large GPU memory, which slowed down training.
1. **Ineffectiveness of fine-tuning**: Efforts to fine-tune the ALIGNN-FF WT10 model on the CHGNet data suffered high initial loss, even worse than the untrained model, indicating significant dataset incompatibility.

The decision to abort testing ALIGNN FF was made after weeks of work due to ongoing technical challenges and resource limitations. See the [PR discussion](https://github.com/janosh/matbench-discovery/pull/47) for further details.

### Fine-tuning

We attempted fine-tuning the [`alignnff_wt10` checkpoint](https://github.com/usnistgov/alignn/blob/461b35fe/alignn/ff/alignnff_wt10/best_model.pt).

The patch [`alignn-ff/alignn-ff-2023.07.05.patch`](alignn-ff/alignn-ff-2023.07.05.patch)
fixes the following issue:

```bash
Traceback (most recent call last):
  File "alignn_relax.py", line 96, in <module>
  File "alignn_relax.py", line 88, in alignn_relax
  File "../alignn/ff/ff.py", line 310, in optimize_atoms
  File "../alignn/lib/python3.9/site-packages/ase/optimize/optimize.py", line 269, in run
  File "../alignn/lib/python3.9/site-packages/ase/optimize/optimize.py", line 156, in run
  File "../alignn/lib/python3.9/site-packages/ase/optimize/optimize.py", line 129, in irun
  File "../alignn/lib/python3.9/site-packages/ase/optimize/optimize.py", line 108, in call_observers
  File "../alignn/lib/python3.9/site-packages/ase/io/trajectory.py", line 132, in write
  File "../alignn/lib/python3.9/site-packages/ase/io/trajectory.py", line 156, in _write_atoms
  File "../alignn/lib/python3.9/site-packages/ase/io/trajectory.py", line 381, in write_atoms
  File "../alignn/lib/python3.9/site-packages/ase/io/ulm.py", line 400, in write
  File "../alignn/lib/python3.9/site-packages/ase/io/ulm.py", line 325, in fill
OSError: [Errno 24] Too many open files
```

The abandoned relaxation and discovery scripts were removed when
[PR #375](https://github.com/janosh/matbench-discovery/pull/375) centralized
discovery execution. This directory retains the metadata and patches as a
historical record of the aborted submission.
