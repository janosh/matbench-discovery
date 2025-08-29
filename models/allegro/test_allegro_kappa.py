"""
Script for generating the predicted kappa-SRME values for the 103 structures
in the PhononDB-PBE dataset, using a NequIP/Allegro model.

Templated from
https://github.com/janosh/matbench-discovery/blob/main/models/mace/calc_kappa_mace_ray_parallelized.py
"""

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1

import contextlib
import json
import os
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from glob import glob
from importlib.metadata import version
from typing import TYPE_CHECKING, Any, Literal

import ase.io
import ase.optimize as opt
import pandas as pd
import torch
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, FrechetCellFilter
from nequip.ase import NequIPCalculator
from pymatgen.core.structure import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from collections.abc import Callable

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

matching_files = glob(f"{compile_path}")
if len(matching_files) == 1:
    compiled_model_file = next(iter(matching_files))
elif os.path.exists(f"{compile_path}"):
    compiled_model_file = f"{compile_path}"
else:
    raise FileNotFoundError(f"No compiled model file was not found at {compile_path}!")


def calc_kappa_for_structure(
    *,  # force keyword-only arguments
    atoms: Atoms,
    displacement_distance: float,
    temperatures: list[float],
    ase_optimizer: str,
    ase_filter: str,
    max_steps: int,
    force_max: float,
    symprec: float,
    enforce_relax_symm: bool,
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
        compile_path (str): File path to compiled model checkpoint
        temperatures (list[float]): Which temperatures to calculate kappa at in Kelvin
        ase_optimizer (str): ASE optimizer to use
        ase_filter (str): ASE filter to use
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
    print(f"Calculating {Structure.from_ase_atoms(atoms).reduced_formula}")
    warnings.filterwarnings("ignore", category=FutureWarning, module="torch")
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using {device=}")

    # Initialize Nequip ASE Calculator from checkpoint
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Trying to use model type names")
        calc = NequIPCalculator.from_compiled_model(
            compile_path=compiled_model_file,
            device="cuda" if torch.cuda.is_available() else "cpu",
        )

    # Ensure arrays are writable
    atoms.arrays = {key: val.copy() for key, val in atoms.arrays.items()}

    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    mat_name = atoms.info["name"]
    info_dict: dict[str, Any] = {
        "name": mat_name,
        "errors": [],
        "error_traceback": [],
    }

    filter_cls: Callable[[Atoms], Atoms] = {
        "frechet": FrechetCellFilter,
        "exp": ExpCellFilter,
    }[ase_filter]
    optimizer_dict = {
        "GPMin": opt.GPMin,
        "GOQN": opt.GoodOldQuasiNewton,
        "BFGSLineSearch": opt.BFGSLineSearch,
        "QuasiNewton": opt.BFGSLineSearch,
        "SciPyFminBFGS": opt.sciopt.SciPyFminBFGS,
        "BFGS": opt.BFGS,
        "LBFGSLineSearch": opt.LBFGSLineSearch,
        "SciPyFminCG": opt.sciopt.SciPyFminCG,
        "FIRE2": opt.fire2.FIRE2,
        "FIRE": opt.fire.FIRE,
        "LBFGS": opt.LBFGS,
    }
    optim_cls: Callable[..., opt.optimize.Optimizer] = optimizer_dict[ase_optimizer]

    # Initialize variables that might be needed in error handling
    relax_dict = {"max_stress": None, "reached_max_steps": False}
    force_results = None

    # Relaxation
    try:
        atoms.calc = calc
        if max_steps > 0:
            if enforce_relax_symm:
                atoms.set_constraint(FixSymmetry(atoms))
                filtered_atoms = filter_cls(atoms, mask=[True] * 3 + [False] * 3)
            else:
                filtered_atoms = filter_cls(atoms)

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
        warnings.warn(f"Failed to relax {mat_name=}, {mat_id=}: {exc!r}", stacklevel=2)
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
        ltc_condition = ignore_imaginary_freqs or not has_imaginary_freqs

        if ltc_condition:
            fc3_set = ltc.calculate_fc3_set(
                ph3, calculator=calc, pbar_kwargs={"position": task_id}
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            warnings.warn(
                f"Imaginary frequencies calculated for {mat_id}, "
                f"skipping FC3 and LTC calculation!",
                stacklevel=2,
            )
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


job_name = f"kappa-103-{ase_optimizer}-dist={displacement_distance}-{fmax=}-{symprec=}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{model_name}/{today}-{job_name}")
os.makedirs(out_dir, exist_ok=True)
timestamp = f"{datetime.now().astimezone():%Y-%m-%d@%H-%M-%S}"
print(f"\nJob {job_name} with {model_name} started {timestamp}")

atoms_list = ase.io.read(DataFiles.phonondb_pbe_103_structures.path, index=":")
atoms_list = sorted(
    atoms_list, key=len
)  # sort by size to get roughly even distribution of comp cost across GPUs
if slurm_array_task_count > 1:
    atoms_list = atoms_list[
        slurm_array_task_id::slurm_array_task_count
    ]  # even distribution of rough comp cost, based on size

# Save run parameters
remote_params = dict(
    # model_name=model_name,
    compile_path=compiled_model_file,
    ase_optimizer=ase_optimizer,
    ase_filter=ase_filter,
    max_steps=max_steps,
    force_max=fmax,
    symprec=symprec,
    enforce_relax_symm=enforce_relax_symm,
    conductivity_broken_symm=conductivity_broken_symm,
    temperatures=temperatures,
    out_dir=out_dir,
    displacement_distance=displacement_distance,
    save_forces=save_forces,
)
run_params = dict(
    **remote_params,
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
        atoms=atoms, **remote_params, task_id=idx
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
