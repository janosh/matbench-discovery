"""AlphaNet thermal conductivity calculation script."""

import json
import os
import traceback
import warnings
from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from typing import TYPE_CHECKING, Any, Literal, cast

import pandas as pd
import torch
from alphanet.config import All_Config
from alphanet.infer.calc import AlphaNetCalculator
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from ase import Atoms

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")

# EDITABLE CONFIG

model_name = "alphanet"
device = "cuda" if torch.cuda.is_available() else "cpu"
dtype = "float64"
config = All_Config().from_json("./mp/mp.json")

calc = AlphaNetCalculator(
    ckpt_path="./mp/mp_0329.ckpt",
    device="cuda",
    precision="64",
    config=config,
)


# Relaxation parameters
ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
ase_filter: Literal["frechet", "exp"] = "frechet"
max_steps = 300
force_max = 1e-4  # Run until the forces are smaller than this in eV/A

# Symmetry parameters
# symmetry precision for enforcing relaxation and conductivity calculation
symprec = 1e-5
# Enforce symmetry with during relaxation if broken
enforce_relax_symm = True
# Conductivity to be calculated if symmetry group changed during relaxation
conductivity_broken_symm = False
prog_bar = True
save_forces = False  # Save force sets to file
temperatures = [300]  # Temperatures to calculate conductivity at in Kelvin
displacement_distance = 0.01  # Displacement distance for phono3py

idx = 1

task_type = "LTC"  # lattice thermal conductivity
job_name = (
    f"{model_name}-phononDB-{task_type}-{ase_optimizer}_force{force_max}_sym{symprec}"
)
module_dir = os.path.dirname(__file__)
out_dir = "./results"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/conductivity_{idx}.json.gz"

timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
struct_data_path = f"../data/part_{idx}.extxyz"
print(f"\nJob {job_name} started {timestamp}")
print(f"Read data from {struct_data_path}")
atoms_list = cast("list[Atoms]", read(struct_data_path, format="extxyz", index=":"))

run_params = {
    "timestamp": timestamp,
    "model_name": model_name,
    "device": device,
    "dtype": dtype,
    "versions": {dep: version(dep) for dep in ("numpy", "torch", "matbench_discovery")},
    "ase_optimizer": ase_optimizer,
    "ase_filter": ase_filter,
    "max_steps": max_steps,
    "force_max": force_max,
    "symprec": symprec,
    "enforce_relax_symm": enforce_relax_symm,
    "conductivity_broken_symm": conductivity_broken_symm,
    "temperatures": temperatures,
    "displacement_distance": displacement_distance,
    "task_type": task_type,
    "job_name": job_name,
    "struct_data_path": os.path.basename(struct_data_path),
    "n_structures": len(atoms_list),
}

with open(f"{out_dir}/run_params.json", "w") as f:
    json.dump(run_params, f, indent=4)

# Set up the relaxation and force set calculation
optim_cls: Callable[..., Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
force_results: dict[str, dict[str, Any]] = {}
kappa_results: dict[str, dict[str, Any]] = {}
tqdm_bar = tqdm(atoms_list, desc="Conductivity calculation: ", disable=not prog_bar)

for atoms in tqdm_bar:
    mat_id = atoms.info.get(Key.mat_id, f"id-{len(kappa_results)}")
    init_info = deepcopy(atoms.info)
    mat_name = atoms.info.get("name", "unknown")

    spg_num = MoyoDataset(MoyoAdapter.from_atoms(atoms)).number
    mat_desc = f"{mat_name}-{spg_num}"

    info_dict = {
        "desc": mat_desc,
        "name": mat_name,
        "initial_space_group_number": spg_num,
        "errors": [],
        "error_traceback": [],
    }

    tqdm_bar.set_postfix_str(mat_desc, refresh=True)

    # Initialize relax_dict to avoid "possibly unbound" errors
    relax_dict = {
        "max_stress": None,
        "reached_max_steps": False,
        "broken_symmetry": False,
    }

    try:
        atoms.calc = calc
        if max_steps > 0:
            if enforce_relax_symm:
                atoms.set_constraint(FixSymmetry(atoms))
                # Use standard mask for no-tilt constraint
                filtered_atoms = FrechetCellFilter(atoms, mask=[True] * 3 + [False] * 3)
            else:
                filtered_atoms = FrechetCellFilter(atoms)

            optimizer = optim_cls(
                filtered_atoms, logfile=f"{out_dir}/relax_{mat_id}.log"
            )
            optimizer.run(fmax=force_max, steps=max_steps)

            reached_max_steps = optimizer.step >= max_steps
            if reached_max_steps:
                print(f"Material {mat_desc=} reached {max_steps=} during relaxation.")

            max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)
            atoms.calc = None
            atoms.constraints = None
            atoms.info = init_info | atoms.info

            # Check if symmetry was broken during relaxation
            relaxed_spg = MoyoDataset(MoyoAdapter.from_atoms(atoms)).number
            broken_symmetry = spg_num != relaxed_spg
            relax_dict = {
                "max_stress": max_stress,
                "reached_max_steps": reached_max_steps,
                "relaxed_space_group_number": relaxed_spg,
                "broken_symmetry": broken_symmetry,
            }

    except Exception as exc:
        warnings.warn(f"Failed to relax {mat_name=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"RelaxError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict
        continue

    # Calculation of force sets
    try:
        # Initialize phono3py with the relaxed structure
        ph3 = ltc.init_phono3py(
            atoms,
            fc2_supercell=atoms.info.get("fc2_supercell", [2, 2, 2]),
            fc3_supercell=atoms.info.get("fc3_supercell", [2, 2, 2]),
            q_point_mesh=atoms.info.get("q_point_mesh", [10, 10, 10]),
            displacement_distance=displacement_distance,
            symprec=symprec,
        )

        # Calculate force constants and frequencies
        ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
            ph3,
            calculator=calc,
            pbar_kwargs={"leave": False, "disable": not prog_bar},
        )

        # Check for imaginary frequencies
        has_imaginary_freqs = check_imaginary_freqs(freqs)
        freqs_dict = {
            Key.has_imag_ph_modes: has_imaginary_freqs,
            Key.ph_freqs: freqs,
        }

        # If conductivity condition is met, calculate fc3
        ltc_condition = not has_imaginary_freqs and (
            not relax_dict["broken_symmetry"] or conductivity_broken_symm
        )

        if ltc_condition:  # Calculate third-order force constants
            print(f"Calculating FC3 for {mat_id}")
            fc3_set = ltc.calculate_fc3_set(
                ph3,
                calculator=calc,
                pbar_kwargs={"leave": False, "disable": not prog_bar},
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            fc3_set = []

        if save_forces:
            force_results[mat_id] = {"fc2_set": fc2_set, "fc3_set": fc3_set}

        if not ltc_condition:
            kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
            warnings.warn(
                f"Material {mat_desc} imaginary frequencies or broken symmetry.",
                stacklevel=2,
            )
            continue

    except Exception as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"ForceConstantError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict
        continue

    try:  # Calculate thermal conductivity
        ph3, kappa_dict, _cond = ltc.calculate_conductivity(
            ph3, temperatures=temperatures
        )
        print(f"Calculated kappa for {mat_id}: {kappa_dict}")
    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        info_dict["errors"].append(f"ConductivityError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
        continue

    kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict

# Save results
df_kappa = pd.DataFrame(kappa_results).T
df_kappa.index.name = Key.mat_id
df_kappa.reset_index().to_json(out_path)
print(f"Saved kappa results to {out_path}")

if save_forces:
    force_out_path = f"{out_dir}/force_sets.json.gz"
    df_force = pd.DataFrame(force_results).T
    df_force.index.name = Key.mat_id
    df_force.reset_index().to_json(force_out_path)
    print(f"Saved force results to {force_out_path}")
