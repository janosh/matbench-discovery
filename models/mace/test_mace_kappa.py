# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "ase",
#   "cuequivariance",
#   "cuequivariance-ops-torch-cu12",
#   "cuequivariance-torch",
#   "mace-torch>=0.3.15",
#   "matbench-discovery[phonons]",
#   "pandas",
#   "pymatviz",
#   "torch",
#   "tqdm",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///

"""This script runs MACE kappa predictions without Ray for consistency."""

import json
import os
import warnings
from datetime import datetime
from importlib.metadata import version
from typing import Any, Final, Literal

import ase.io
import pandas as pd
import torch
from mace.calculators import mace_mp
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.enums import DataFiles
from matbench_discovery.phonons import KappaCalcParams, calc_kappa_for_structure

module_dir = os.path.dirname(__file__)

# Relaxation parameters
ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
max_steps = 300
force_max = 1e-4  # Run until the forces are smaller than this in eV/A

# Symmetry parameters
symprec = 1e-5  # symmetry precision for enforcing relaxation and conductivity calcs
enforce_relax_symm = True  # Enforce symmetry with during relaxation if broken
# Conductivity to be calculated if symmetry group changed during relaxation
conductivity_broken_symm = False
save_forces = True  # Save force sets to file
temperatures: list[float] = [300]

model_name = os.getenv("MODEL_NAME", "mace_omat_0")
release_base = "https://github.com/ACEsuit/mace-foundations/releases/download"
model_configs: Final[dict[str, dict[str, str]]] = {
    "mace_omat_0": {
        "checkpoint": f"{release_base}/mace_omat_0/mace-omat-0-medium.model",
    },
    "mace_matpes_0": {
        "checkpoint": f"{release_base}/mace_matpes_0/MACE-matpes-pbe-omat-ft.model",
    },
    "mace_mh_1": {
        "checkpoint": f"{release_base}/mace_mh_1/mace-mh-1.model",
        "head": "omat_pbe",
    },
}
model_config = model_configs[model_name]
checkpoint = model_config["checkpoint"]

print(f"Loading {model_name}...")
warnings.filterwarnings("ignore", category=FutureWarning, module="torch")
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device=}")
mace_kwargs: dict[str, Any] = {}
if "head" in model_config:
    mace_kwargs["head"] = model_config["head"]
mace_calc = mace_mp(
    model=checkpoint, device=device, enable_cueq=device == "cuda", **mace_kwargs
)

displacement_distance = 0.01
job_name = (
    f"kappa-103-{ase_optimizer}-dist={displacement_distance}-"
    f"fmax={force_max}-{symprec=}"
)
out_dir = f"{module_dir}/{model_name}/{today}-{job_name}"
os.makedirs(out_dir, exist_ok=True)

timestamp = f"{datetime.now().astimezone():%Y-%m-%d@%H-%M-%S}"
print(f"\nJob {job_name} with {model_name} started {timestamp}")

atoms_list = ase.io.read(DataFiles.phonondb_pbe_103_structures.path, index=":")

# Save run parameters
kappa_params: KappaCalcParams = {
    "ase_optimizer": ase_optimizer,
    "max_steps": max_steps,
    "force_max": force_max,
    "symprec": symprec,
    "enforce_relax_symm": enforce_relax_symm,
    "conductivity_broken_symm": conductivity_broken_symm,
    "temperatures": temperatures,
    "out_dir": out_dir,
    "displacement_distance": displacement_distance,
    "save_forces": save_forces,
}
run_params = dict(
    **kappa_params,
    model_name=model_name,
    checkpoint=checkpoint,
    mace_kwargs=mace_kwargs,
    n_structures=len(atoms_list),
    struct_data_path=DataFiles.phonondb_pbe_103_structures.path,
    versions={dep: version(dep) for dep in ("numpy", "torch", "mace-torch")},
)

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Process results as they complete
kappa_results: dict[str, dict[str, Any]] = {}
force_results: dict[str, dict[str, Any]] = {}

for idx, atoms in enumerate(
    tqdm(
        atoms_list,
        desc=f"{model_name} kappa",
        total=len(atoms_list),
        mininterval=5,
        unit="struct",
    )
):
    mat_id, result_dict, force_dict = calc_kappa_for_structure(
        atoms=atoms,
        calculator=mace_calc,
        **kappa_params,
        task_id=idx,
    )
    kappa_results[mat_id] = result_dict
    if force_dict is not None:
        force_results[mat_id] = force_dict

    # Save intermediate results
    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    if Key.mat_id not in df_kappa:
        df_kappa = df_kappa.reset_index()
    df_kappa.to_json(f"{out_dir}/kappa.json.gz")

    if save_forces:
        df_force = pd.DataFrame(force_results).T
        df_force = df_force.drop(columns=Key.mat_id, errors="ignore")
        df_force = pd.concat([df_kappa.set_index(Key.mat_id), df_force], axis=1)
        df_force.index.name = Key.mat_id
        df_force.reset_index().to_json(f"{out_dir}/force-sets.json.gz")

print(f"\nResults saved to {out_dir!r}")
