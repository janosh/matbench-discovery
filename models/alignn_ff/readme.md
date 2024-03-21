# ALIGNN-FF (submission aborted)

The [ALIGNN FF model submission](https://github.com/janosh/matbench-discovery/pull/47) intended to get a complete set of formation energy predictions for the WBM test set post-ALIGNN-FF structure relaxation (i.e. the WBM IS2RE task).

This effort was aborted for the following reasons:

1. **Incompatibility issues**: ALIGNN-FF was pre-trained on the JARVIS data, which among other differences uses the OptB88vdW functional and is incompatible with the WBM test set generated using Materials Project workflows.
1. **Training difficulties**: ALIGNN-FF proved to be very resource-hungry. [12 GB of MPtrj training data](https://figshare.com/articles/dataset/23713842) turned into 600 GB of ALIGNN graph data. This forces small batch size even on nodes with large GPU memory, which slowed down training.
1. **Ineffectiveness of fine-tuning**: Efforts to fine-tune the ALIGNN-FF WT10 model on the CHGNet data suffered high initial loss, even worse than the untrained model, indicating significant dataset incompatibility.

The decision to abort testing ALIGNN FF was made after weeks of work due to ongoing technical challenges and resource limitations. See the [PR discussion](https://github.com/janosh/matbench-discovery/pull/47) for further details.

## Fine-tuning

We attempted fine-tuning the [`alignnff_wt10` checkpoint](https://github.com/usnistgov/alignn/blob/461b35fe6e5ed7ade7cbf9b345773e941371ecfc/alignn/ff/alignnff_wt10/best_model.pt).

The patch `alignn-ff-2023.07.05.patch` fixes the following issue:

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

## Scripts

1. `alignn_ff_relax.py`: Relax WBM test set structures. Set the variable `n_splits` to the number of GPU compute nodes. On each compute node, set the environment variable `TASK_ID` to a value in the range 1-`n_splits`. Set the variable `n_processes_per_task` to the number of processes on a single node. For 48 CPU cores with 4 GPUs a good setting is to use 10 processes.
2. `test_alignn_ff.py`: Read the relaxed structures from `alignn_ff_relax.py` and make formation energy predictions. Set the variable `n_splits` accordingly.
