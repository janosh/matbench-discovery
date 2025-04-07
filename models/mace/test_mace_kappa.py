"""This script runs MACE predictions in parallel using Ray.
Can be run both locally and on a cluster with automatic resource scaling.
"""

import json
import os
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from typing import TYPE_CHECKING, Any, Final, Literal

import ase.io
import pandas as pd
import ray
import torch
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from mace.calculators import mace_mp
from pymatviz.enums import Key
from tqdm import tqdm

import matbench_discovery
from matbench_discovery import today
from matbench_discovery.enums import DataFiles
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from ase.optimize.optimize import Optimizer

module_dir = os.path.dirname(__file__)


@ray.remote(
    num_cpus=1,
    # num_gpus=1,
)
def calc_kappa_for_structure(
    *,  # force keyword-only arguments
    atoms: Atoms,
    displacement_distance: float,
    checkpoint: str,
    temperatures: list[float],
    ase_optimizer: str,
    max_steps: int,
    force_max: float,
    symprec: float,
    enforce_relax_symm: bool,
    conductivity_broken_symm: bool,
    save_forces: bool,
    out_dir: str,
    task_id: int,
    device: str | None = None,
) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    """Predict ML kappa for single structure with ray.

    Args:
        atoms (Atoms): ASE Atoms object with fc2_supercell, fc3_supercell,
            q_point_mesh keys in its info dict.
        displacement_distance (float): Displacement distance for phono3py
        checkpoint (str): File path or download URL to model checkpoint
        temperatures (list[float]): Which temperatures to calculate kappa at in Kelvin
        ase_optimizer (str): ASE optimizer to use
        max_steps (int): Maximum number of optimization steps
        force_max (float): Maximum force tolerance
        symprec (float): Symmetry precision
        enforce_relax_symm (bool): Whether to enforce symmetry during relaxation
        conductivity_broken_symm (bool): Whether to calculate conductivity if
            symmetry broken
        save_forces (bool): Whether to save force sets
        out_dir (str): Output directory
        task_id (int): Task ID for logging

    Returns:
        tuple[str, dict[str, Any], dict[str, Any] | None]:
            material ID, results dict, force results dict
    """
    warnings.filterwarnings("ignore", category=FutureWarning, module="torch")
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using {device=}")
    calc = mace_mp(model=checkpoint, device=device, enable_cueq=device == "cuda")

    # Create a deep copy of the atoms object to avoid ray read-only issues
    atoms = atoms.copy()
    # Ensure arrays are writable
    atoms.arrays = {key: val.copy() for key, val in atoms.arrays.items()}

    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    formula = atoms.get_chemical_formula()
    info_dict: dict[str, Any] = {
        Key.formula: formula,
        "errors": [],
        "error_traceback": [],
    }
    optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

    # Initialize variables that might be needed in error handling
    relax_dict = {"max_stress": None, "reached_max_steps": False}
    force_results = None

    # Relaxation
    try:
        atoms.calc = calc
        if max_steps > 0:
            if enforce_relax_symm:
                atoms.set_constraint(FixSymmetry(atoms))
                filtered_atoms = FrechetCellFilter(atoms, mask=[True] * 3 + [False] * 3)
            else:
                filtered_atoms = FrechetCellFilter(atoms)

            os.makedirs(relax_dir := f"{out_dir}/relaxations", exist_ok=True)
            optimizer = optim_cls(filtered_atoms, logfile=f"{relax_dir}/{task_id}.log")
            optimizer.run(fmax=force_max, steps=max_steps)

            reached_max_steps = optimizer.step == max_steps
            if reached_max_steps:
                print(f"Material {mat_id=} reached {max_steps=} during relaxation")

            # maximum residual stress component in for xx,yy,zz and xy,yz,xz
            # components separately result is a array of 2 elements
            max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

            atoms.calc = None
            atoms.constraints = None
            atoms.info = init_info | atoms.info

            relax_dict = {
                "max_stress": max_stress,
                "reached_max_steps": reached_max_steps,
            }

    except Exception as exc:
        warnings.warn(f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"] += [f"RelaxError: {exc!r}"]
        info_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict, None

    # Calculation of force sets
    try:
        ph3 = ltc.init_phono3py(
            atoms,
            fc2_supercell=atoms.info["fc2_supercell"],
            fc3_supercell=atoms.info["fc3_supercell"],
            q_point_mesh=atoms.info["q_point_mesh"],
            displacement_distance=displacement_distance,
            symprec=symprec,
        )

        ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
            ph3, calculator=calc, pbar_kwargs={"disable": True}
        )

        has_imaginary_freqs = check_imaginary_freqs(freqs)
        freqs_dict = {Key.has_imag_ph_modes: has_imaginary_freqs, Key.ph_freqs: freqs}

        # if conductivity condition is met, calculate fc3
        ltc_condition = not has_imaginary_freqs and (
            not relax_dict["broken_symmetry"] or conductivity_broken_symm
        )

        if ltc_condition:
            fc3_set = ltc.calculate_fc3_set(
                ph3, calculator=calc, pbar_kwargs={"position": task_id}
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            fc3_set = []

        force_results = (
            {"fc2_set": fc2_set, "fc3_set": fc3_set} if save_forces else None
        )

        if not ltc_condition:
            return mat_id, info_dict | relax_dict | freqs_dict, force_results

    except Exception as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"] += [f"ForceConstantError: {exc!r}"]
        info_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict, force_results

    # Calculation of conductivity
    try:
        ph3, kappa_dict, _cond = ltc.calculate_conductivity(
            ph3, temperatures=temperatures
        )
        return mat_id, info_dict | relax_dict | freqs_dict | kappa_dict, force_results

    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        info_dict["errors"] += [f"ConductivityError: {exc!r}"]
        info_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict | freqs_dict, force_results


# Relaxation parameters
ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
max_steps = 300
fmax = 1e-4  # Run until the forces are smaller than this in eV/A

# Symmetry parameters
symprec = 1e-5  # symmetry precision for enforcing relaxation and conductivity calcs
enforce_relax_symm = True  # Enforce symmetry with during relaxation if broken
# Conductivity to be calculated if symmetry group changed during relaxation
conductivity_broken_symm = False
save_forces = True  # Save force sets to file
temperatures: list[float] = [300]

# Initialize Ray (point this at the head node of the cluster)
ray_ip: Final[str] = {
    "lambda-staging-with-ray-2.40": "100.82.154.22",
}.get(os.getenv("RAY_IP_KEY", "")) or "100.82.154.22"
ray_address = os.getenv("RAY_ADDRESS", f"ray://{ray_ip}:10001")
print(f"{ray_address=}")

# Ray initialization
if ray_address:
    # Connect to existing Ray cluster
    ray.init(
        address=ray_address,
        runtime_env={
            "py_modules": [matbench_discovery],
            # "working_dir": PKG_DIR,  # add matbench-discovery root to PYTHONPATH
            "uv": [
                # "mace-torch",
                "git+https://github.com/janosh/mace",
                # "cuequivariance-ops-torch-cu12",
                # "cuequivariance-torch",
                # "cuequivariance",
                "phono3py",
                "ase",
                "moyopy",
                "pymatviz",
            ],
        },
    )
else:
    # Start Ray locally with optimized settings for M3 Max
    ray.init(num_cpus=8, num_gpus=1)

print(f"\nConnected to Ray cluster: {ray.cluster_resources()}")
ray_resources = ray.available_resources()
ray_mem = ray_resources.get("memory", 0) / 1e9
print(f"Available memory: {ray_mem:.1f} GB")
obj_store_mem = ray_resources.get("object_store_memory", 0)
print(f"Object store memory: {obj_store_mem / 1e9:.1f} GB")

model_name = "mace-omat-0-medium"
checkpoint = f"https://github.com/ACEsuit/mace-foundations/releases/download/mace_omat_0/{model_name}.model"

displacement_distance = 0.01
job_name = (
    f"{today}-kappa-103-{ase_optimizer}-dist={displacement_distance}-{fmax=}-{symprec=}"
)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{model_name}/{job_name}")
os.makedirs(out_dir, exist_ok=True)

timestamp = f"{datetime.now().astimezone():%Y-%m-%d@%H-%M-%S}"
print(f"\nJob {job_name} with {model_name} started {timestamp}")

atoms_list = ase.io.read(DataFiles.phonondb_pbe_103_structures.path, index=":")
# Save run parameters
remote_params = dict(
    model_name=model_name,
    checkpoint=checkpoint,
    ase_optimizer=ase_optimizer,
    max_steps=max_steps,
    force_max=fmax,
    symprec=symprec,
    enforce_relax_symm=enforce_relax_symm,
    conductivity_broken_symm=conductivity_broken_symm,
    temperatures=temperatures,
    out_dir=out_dir,
)
run_params = dict(
    **remote_params,
    n_structures=len(atoms_list),
    struct_data_path=DataFiles.phonondb_pbe_103_structures.path,
    versions={dep: version(dep) for dep in ("numpy", "torch", "ray")},
)

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Process structures in parallel
futures = [
    calc_kappa_for_structure.remote(atoms=atoms, **remote_params, task_id=idx)
    for idx, atoms in enumerate(atoms_list[:])
]

# Process results as they complete
kappa_results: dict[str, dict[str, Any]] = {}
force_results: dict[str, dict[str, Any]] = {}

for future in tqdm(futures, desc=f"Predicting kappa with {model_name}"):
    mat_id, result_dict, force_dict = ray.get(future)
    kappa_results[mat_id] = result_dict
    if force_dict is not None:
        force_results[mat_id] = force_dict

    # Save intermediate results
    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    df_kappa.reset_index().to_json(f"{out_dir}/kappa.json.gz")

    if save_forces:
        df_force = pd.DataFrame(force_results).T
        df_force = pd.concat([df_kappa, df_force], axis=1)
        df_force.index.name = Key.mat_id
        df_force.reset_index().to_json(f"{out_dir}/force-sets.json.gz")

print(f"\nResults saved to {out_dir!r}")
