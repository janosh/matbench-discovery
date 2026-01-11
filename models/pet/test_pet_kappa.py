"""
Script for generating the predicted kappa-SRME values for the 103 structures in the
PhononDB-PBE dataset, using a PET model.

Templated from https://github.com/janosh/matbench-discovery/blob/main/models/nequip/test_nequip_kappa.py
"""

import json
import os
import traceback
import warnings
from datetime import datetime
from importlib.metadata import version
from typing import Any, Literal

import ase.io
import pandas as pd
import torch
from calc_kappa import calc_kappa_for_structure
from metatomic.torch import load_atomistic_model
from metatomic.torch.ase_calculator import MetatomicCalculator, SymmetrizedCalculator
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles
from matbench_discovery.metrics.phonons import calc_kappa_metrics_from_dfs
from matbench_discovery.phonons import KappaCalcParams

# Model configuration
module_dir = os.path.dirname(__file__)
model_name = "pet"
model_variant = "oam-xl-v1.0.0"  # get it with `mtt export https://huggingface.co/lab-cosmo/upet/resolve/main/models/pet-oam-xl-v1.0.0.ckpt`
precision = "float64"
device = "cuda" if torch.cuda.is_available() else "cpu"
dtype = torch.float64 if precision == "float64" else torch.float32
model = load_atomistic_model(f"{model_name}-{model_variant}.pt")
model.capabilities().dtype = precision
model = model.to(dtype=dtype, device=device)
calc = MetatomicCalculator(model, device=device, non_conservative=False)
calc = SymmetrizedCalculator(calc, batch_size=1, include_inversion=False)
batch_size = 1

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
    versions={dep: version(dep) for dep in ("numpy", "torch", "metatomic")},
)

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Process results as they complete
kappa_results: dict[str, dict[str, Any]] = {}
force_results: dict[str, dict[str, Any]] = {}

for idx, atoms in enumerate(tqdm(atoms_list, desc="Calculating kappa...")):
    mat_id, result_dict, force_dict = calc_kappa_for_structure(
        atoms=atoms,
        calculator=calc,
        batch_size=batch_size,
        is_plusminus=True,
        ignore_imaginary_freqs=ignore_imaginary_freqs,
        formula_getter=lambda a: a.info.get("name", a.get_chemical_formula()),
        **kappa_params,
        task_id=idx,
    )
    kappa_results[mat_id] = result_dict
    if force_dict is not None:
        force_results[mat_id] = force_dict

    # Save intermediate results
    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    df_kappa.reset_index(drop=True).to_json(
        f"{out_dir}/{slurm_array_task_id}_kappa.json.gz"
    )
    df_kappa.to_json(f"{out_dir}/{slurm_array_task_id}_kappa.json.gz")

    if save_forces:
        df_force = pd.DataFrame(force_results).T
        df_force = pd.concat([df_kappa, df_force], axis=1)
        df_force.index.name = Key.mat_id
        df_force.reset_index(drop=True).to_json(
            f"{out_dir}/{slurm_array_task_id}_force-sets.json.gz"
        )

print(f"\nResults saved to {out_dir!r}")

try:
    print("Computing metrics against reference data...")
    df_dft = pd.read_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path).set_index(
        Key.mat_id
    )
    if ignore_imaginary_freqs:
        # WARNING: setting has_imag_ph_modes to False to compute the metrics anyway
        df_kappa["has_imag_ph_modes"] = False
    df_ml_metrics = calc_kappa_metrics_from_dfs(df_kappa, df_dft)
    # Compute and print summary metrics
    kappa_sre = df_ml_metrics[Key.sre].mean()
    kappa_srme = df_ml_metrics[Key.srme].mean()
    print(f"{kappa_sre=:.4f}")
    print(f"{kappa_srme=:.4f}")

    df_ml_metrics.to_json(f"{out_dir}/metrics.json.gz")
    print(f"Saved metrics to {out_dir}/metrics.json.gz")
except Exception as exc:
    warnings.warn(f"Failed to calculate metrics: {exc!r}", stacklevel=2)
    traceback.print_exc()
