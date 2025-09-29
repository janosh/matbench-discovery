"""
Script for testing predictions of a trained NequIP/Allegro model on the WBM
test dataset (used as the matbench-discovery test set for models trained on
MPTrj).
Copied from the 7net script here:
https://github.com/janosh/matbench-discovery/blob/main/models/sevennet/test_sevennet.py
Then refactored for NequIP/Allegro and SLURM on the Frontier HPC (allowing
parallelization over many GPUs/nodes).
"""

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1

import contextlib
import os
import warnings
from glob import glob
from typing import Any, Literal

import ase.optimize
import ase.optimize.sciopt
import numpy as np
import pandas as pd
import torch
from ase.filters import ExpCellFilter, Filter, FrechetCellFilter
from ase.optimize import Optimizer
from nequip.ase import NequIPCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import timestamp
from matbench_discovery.data import DataFiles, as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import Task

with contextlib.suppress(ImportError):
    # OpenEquivariance/CuEquivariance libraries need to be loaded to allow their use
    # in ASE calculators, if model was compiled with these accelerations
    # (see NequIP/Allegro docs), so here we try to import them in case models were
    # compiled with these settings
    pass


# %% this config is editable
compile_path = "*.nequip.pt2"
smoke_test = False  # True
model_name = "allegro-0"
task_type = Task.IS2RE
ase_optimizer = "GOQN"  # faster than "FIRE" from tests, gives the same results;
# see SI of https:/doi.org/10.1088/2515-7655/ade916
ase_filter: Literal["frechet", "exp"] = "frechet"  # recommended filter

max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A

slurm_nodes = int(os.getenv("SLURM_NNODES", "1"))
slurm_tasks_per_node = int(os.getenv("SLURM_NTASKS_PER_NODE", "1"))
slurm_array_task_count = int(os.getenv("NGPUS", slurm_nodes * slurm_tasks_per_node))
slurm_array_task_id = int(
    os.getenv(
        "TASK_ID", os.getenv("SLURM_ARRAY_TASK_ID", os.getenv("SLURM_PROCID", "0"))
    )
)
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", os.getenv("SLURM_JOBID", "debug"))

# Note that we can also manually override some slurm IDs here if we need to rerun
# just a single subset that failed on a previous eval run, for any reason, setting
# job_id to 0, task_id to the failed task, and task_count to match
# whatever the previous task count was (to ensure the same data splitting):
# slurm_array_job_id = 0
# slurm_array_task_id = 104
# slurm_array_task_count = 128

matching_files = glob(f"{compile_path}")
if len(matching_files) == 1:
    compiled_model_file = next(iter(matching_files))
elif os.path.exists(f"{compile_path}"):
    compiled_model_file = f"{compile_path}"
else:
    raise FileNotFoundError(f"No compiled model file was not found at {compile_path}!")

os.makedirs(out_dir := "./results", exist_ok=True)
out_path = f"{out_dir}/{model_name}-{slurm_array_task_id:>03}.json.gz"
job_name = f"{model_name}-wbm-{task_type}-{slurm_array_task_id:>03}"

data_path = {Task.IS2RE: DataFiles.wbm_initial_atoms.path}[
    task_type
]  # automatically downloaded if not already present in cache
print(f"\nJob {job_name!r} running {timestamp}", flush=True)
print(f"{data_path=}", flush=True)
print(f"{slurm_array_task_id} of {slurm_array_task_count}")

# Initialize ASE Calculator from checkpoint
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Trying to use model type names")
    calculator = NequIPCalculator.from_compiled_model(
        compile_path=compiled_model_file,
        device="cuda" if torch.cuda.is_available() else "cpu",
    )


# %%
print(f"Read data from {data_path}")
atoms_list = ase_atoms_from_zip(data_path)
atoms_list = sorted(
    atoms_list, key=len
)  # sort by size to get roughly even distribution of comp cost across GPUs

if slurm_array_job_id == "debug":  # if running a quick smoke test
    if smoke_test:
        atoms_list = atoms_list[:128]
    else:
        pass
elif slurm_array_task_count > 1:
    atoms_list = atoms_list[
        slurm_array_task_id::slurm_array_task_count
    ]  # even distribution of rough comp cost, based on size

relax_results: dict[str, dict[str, Any]] = {}

filter_cls: type[Filter] = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}[ase_filter]
optimizer_dict = {
    "GPMin": ase.optimize.GPMin,
    "GOQN": ase.optimize.GoodOldQuasiNewton,
    "BFGSLineSearch": ase.optimize.BFGSLineSearch,
    "QuasiNewton": ase.optimize.BFGSLineSearch,
    "SciPyFminBFGS": ase.optimize.sciopt.SciPyFminBFGS,
    "BFGS": ase.optimize.BFGS,
    "LBFGSLineSearch": ase.optimize.LBFGSLineSearch,
    "SciPyFminCG": ase.optimize.sciopt.SciPyFminCG,
    "FIRE2": ase.optimize.FIRE2,
    "FIRE": ase.optimize.FIRE,
    "LBFGS": ase.optimize.LBFGS,
}
optim_cls: type[Optimizer] = optimizer_dict[ase_optimizer]


# %%
for atoms in tqdm(atoms_list, desc="Relaxing"):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = calculator
        if max_steps > 0:
            atoms = filter_cls(atoms)
            with optim_cls(atoms, logfile="/dev/null") as optimizer:
                for _ in optimizer.irun(fmax=force_max, steps=max_steps):
                    forces = atoms.get_forces()
                    if np.max(np.linalg.norm(forces, axis=1)) > 1e6:
                        raise RuntimeError(
                            "Forces are exorbitant, exploding relaxation!"
                        )

        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by filter_cls, so extract with getattr
        relaxed_struct = AseAtomsAdaptor.get_structure(getattr(atoms, "atoms", atoms))
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue

df_out = pd.DataFrame(relax_results).T.add_prefix("allegro_")
df_out.index.name = Key.mat_id


# %%
if not smoke_test:
    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)
