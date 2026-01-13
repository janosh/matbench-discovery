"""
Templated from test_nequip_kappa.py and test_orb_kappa.py.
"""

# uses commits matbench-discovery f0e54b7, pymatviz v0.17.3, phono3py 3.22.0

import json
import os
import warnings
from datetime import datetime
from importlib.metadata import version
from typing import Any, Literal

import ase.io
import pandas as pd
import torch
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles, Model
from matbench_discovery.phonons import KappaCalcParams, calc_kappa_for_structure
from matbench_discovery.metrics.phonons import calc_kappa_metrics_from_dfs, write_metrics_to_yaml
from tace.interface.ase import TACEAseCalc

dtype = "float64" 
device = "cuda" if torch.cuda.is_available() else "cpu"
model_name = "TACE-v1-OAM-M"

try:
    from tace.foundations import tace_foundations
    model_path = tace_foundations[model_name]
except Exception as e:
    raise RuntimeError(
        f"Failed to load {model_name}.\n"
        f"Please manual download the model from:\n"
        f"https://huggingface.co/xvzemin/tace-foundations/"
        f"resolve/main/{model_name}.pt\n"
        f"and put the model into ~/.cache/tace/{model_name}.pt"
    ) from e

eval_model = Model.tace_v1_oam_m
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

# Note that we can also manually override some slurm IDs here if we need to rerun just a
# single subset that failed on a previous eval run, for any reason, setting job_id to 0,
# task_id to the failed task, and task_count to match
# whatever the previous task count was (to ensure the same data splitting):
# slurm_array_job_id = 0
# slurm_array_task_id = 104
# slurm_array_task_count = 128

print("Loading TACE model...")
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Trying to use model type names")
    tace_calc = TACEAseCalc(
        model_path,
        use_ema=True,
        device=device,
        dtype=dtype,
    )

job_name = f"kappa-103-{ase_optimizer}-dist={displacement_distance}-{fmax=}-{symprec=}"
out_dir = os.getenv(
    "SBATCH_OUTPUT",
    f"{os.path.dirname(__file__)}/{model_name}/{today}-{job_name}",
)
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
kappa_params: KappaCalcParams = {
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
    versions={dep: version(dep) for dep in ("numpy", "torch", "tace")},
)

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Process results as they complete
kappa_results: dict[str, dict[str, Any]] = {}
force_results: dict[str, dict[str, Any]] = {}

for idx, atoms in enumerate(tqdm(atoms_list, desc="Calculating kappa...")):
    mat_id, result_dict, force_dict = calc_kappa_for_structure(
        atoms=atoms,
        calculator=tace_calc,
        ignore_imaginary_freqs=ignore_imaginary_freqs,
        formula_getter=lambda a: a.info.get("name", a.get_chemical_formula()),
        **kappa_params,
        task_id=idx,
    ) # mat_id like mp-7631
    kappa_results[mat_id] = result_dict
    if force_dict is not None:
        force_results[mat_id] = force_dict

# Save intermediate results
df_kappa = pd.DataFrame(kappa_results).T
df_kappa.index.name = Key.mat_id
df_kappa.reset_index().to_json(f"{out_dir}/{slurm_array_task_id}_kappa.json.gz")

if save_forces:
    df_force = pd.DataFrame(force_results).T
    df_force = pd.concat([pd.DataFrame(kappa_results).T, df_force], axis=1)
    df_force.index.name = Key.mat_id
    df_force.reset_index().to_json(
        f"{out_dir}/{slurm_array_task_id}_force-sets.json.gz"
    )
