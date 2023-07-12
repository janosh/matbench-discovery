# ALIGNN-FF formation energy predictions on WBM test set after ML relaxation

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

To reproduce the ALIGNN relaxed predictions, run the following scripts:

1. `alignn_relax.py`: Set the variable `n_splits` to the number of GPU compute nodes. On each compute node, set the environment variable `TASK_ID` to a value in the range 1-`n_splits`. Set the variable `n_processes_per_task` to the number of processes on a single node. For 48 CPU cores with 4 GPUs a good setting is to use 10 processes.
2. `test_alignn_relaxed.py`: Read the relaxed structures and compute predictions. Set the variable `n_splits` accordingly.
