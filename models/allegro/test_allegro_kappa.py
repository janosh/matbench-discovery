"""
Script for generating the predicted kappa-SRME values for the 103 structures
in the PhononDB-PBE dataset, using a NequIP/Allegro model.

Templated from
https://github.com/janosh/matbench-discovery/blob/main/models/mace/calc_kappa_mace_ray_parallelized.py
"""

# uses commits matbench-discovery 012ccfe, k_srme commit 0269a946, pymatviz v0.15.1

import contextlib
import json
import os
import warnings
from datetime import datetime
from glob import glob
from importlib.metadata import version
from typing import Any, Literal

import ase.io
import pandas as pd
import torch
from nequip.ase import NequIPCalculator
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles
from matbench_discovery.phonons import calc_kappa_for_structure

with contextlib.suppress(ImportError):
    # OpenEquivariance/CuEquivariance libraries need to be loaded to allow their use
    # in ASE calculators, if model was compiled with these accelerations
    # (see NequIP/Allegro docs), so here we try to import them in case models were
    # compiled with these settings
    pass

module_dir = os.path.dirname(__file__)
compile_path = "*.nequip.pt2"
model_name = "allegro-0"

# Relaxation parameters
ase_optimizer = "FIRE"
ase_filter: Literal["frechet", "exp"] = "frechet"  # recommended filter
max_steps = 300
fmax = 1e-4  # Run until the forces are smaller than this in eV/A

# Symmetry parameters
symprec = 1e-5  # symmetry precision for enforcing relaxation and conductivity calcs
enforce_relax_symm = True  # Enforce symmetry with during relaxation if broken
# Conductivity to be calculated if symmetry group changed during relaxation
conductivity_broken_symm = False
save_forces = True  # Save force sets to file
temperatures: list[float] = [300]
displacement_distance = 0.03
ignore_imaginary_freqs = True

# Task splitting:
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

matching_files = glob(compile_path)
if len(matching_files) == 1:
    compiled_model_file = next(iter(matching_files))
elif os.path.isfile(compile_path):
    compiled_model_file = compile_path
else:
    raise FileNotFoundError(f"Compiled model file not found at {compile_path}!")

# Initialize calculator once for all structures
print("Loading Allegro model...")
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Trying to use model type names")
    allegro_calc = NequIPCalculator.from_compiled_model(
        compile_path=compiled_model_file,
        device="cuda" if torch.cuda.is_available() else "cpu",
    )


job_name = f"kappa-103-{ase_optimizer}-dist={displacement_distance}-{fmax=}-{symprec=}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{model_name}/{today}-{job_name}")
os.makedirs(out_dir, exist_ok=True)
timestamp = f"{datetime.now().astimezone():%Y-%m-%d@%H-%M-%S}"
print(f"\nJob {job_name} with {model_name} started {timestamp}")

atoms_list = ase.io.read(DataFiles.phonondb_pbe_103_structures.path, index=":")
# sort by size to get roughly even distribution of comp cost across GPUs
atoms_list = sorted(atoms_list, key=len)
if slurm_array_task_count > 1:
    # even distribution of rough comp cost, based on size
    atoms_list = atoms_list[slurm_array_task_id::slurm_array_task_count]

# Save run parameters
kappa_params = {
    "ase_optimizer": ase_optimizer,
    "ase_filter": ase_filter,
    "max_steps": max_steps,
    "force_max": fmax,
    "symprec": symprec,
    "enforce_relax_symm": enforce_relax_symm,
    "temperatures": temperatures,
    "out_dir": out_dir,
    "displacement_distance": displacement_distance,
    "save_forces": save_forces,
}
run_params = dict(
    **kappa_params,
    n_structures=len(atoms_list),
    struct_data_path=DataFiles.phonondb_pbe_103_structures.path,
    versions={dep: version(dep) for dep in ("numpy", "torch", "nequip")},
)

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Process results as they complete
kappa_results: dict[str, dict[str, Any]] = {}
force_results: dict[str, dict[str, Any]] = {}

for idx, atoms in enumerate(tqdm(atoms_list, desc="Calculating kappa...")):
    mat_id, result_dict, force_dict = calc_kappa_for_structure(
        atoms=atoms,
        calculator=allegro_calc,
        ignore_imaginary_freqs=ignore_imaginary_freqs,
        **kappa_params,  # type: ignore[arg-type]
        task_id=idx,
    )
    kappa_results[mat_id] = result_dict
    if force_dict is not None:
        force_results[mat_id] = force_dict

    # Save intermediate results
    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    df_kappa.reset_index().to_json(f"{out_dir}/{slurm_array_task_id}_kappa.json.gz")

    if save_forces:
        df_force = pd.DataFrame(force_results).T
        df_force = pd.concat([df_kappa, df_force], axis=1)
        df_force.index.name = Key.mat_id
        df_force.reset_index().to_json(
            f"{out_dir}/{slurm_array_task_id}_force-sets.json.gz"
        )

print(f"\nResults saved to {out_dir!r}")
