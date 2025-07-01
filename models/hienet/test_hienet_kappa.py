# /// script
# dependencies = [
#   "torch>=2.1.2",
#   "torch-geometric>=2.6.1",
#   "numpy>=1.26.4",
#   "ase>=3.25.0",
#   "braceexpand>=0.1.7",
#   "e3nn>=0.5.6",
#   "pymatviz>=0.16.0",
#   "pyyaml>=6.0.1",
#   "torch-scatter>=2.1.2",
#   "scikit-learn>=1.7.0",
#   "pymatgen>=2025.6.14",
#   "wandb>=0.20.1",
#   "torch-ema>=0.3",
#   "hienet>=1.0.1",
#   "matbench-discovery>=1.3.1",
#   "phono3py>=3.17.0"
# ]
# ///

import json
import os
import sys
import traceback
import warnings
import argparse
from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from typing import Any, Literal

import pandas as pd
import torch
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from tqdm import tqdm

# Parse command line arguments
parser = argparse.ArgumentParser(description="Thermal conductivity calculation script")
parser.add_argument("--gpu", type=str, default="0", help="GPU ID to use")
parser.add_argument("--left", type=int, default=0, help="Start index for atoms list")
parser.add_argument("--right", type=int, default=None, help="End index for atoms list")
args = parser.parse_args()

from hienet.hienet_calculator import HIENetCalculator

from matbench_discovery import today
from matbench_discovery.enums import DataFiles
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")

# EDITABLE CONFIG
model_name = "HIENet-V3.pth"

os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu
device = "cuda" if torch.cuda.is_available() else "cpu"
model = f"./{model_name}"
calc = HIENetCalculator(model=model, device="cuda")

# Relaxation parameters
ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
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
task_type = "LTC"  # lattice thermal conductivity
job_name = (
    f"{model_name}-phononDB-{task_type}-{ase_optimizer}_force{force_max}_sym{symprec}"
)
module_dir = os.path.dirname(__file__)
out_dir = "./results"
os.makedirs(out_dir, exist_ok=True)
out_path = (
    f"{out_dir}/{today}-HIENet-XL-kappa-103-{ase_optimizer}-dist={displacement_distance}-"
    f"fmax={force_max}-{symprec=}-gpu{args.gpu}-{args.left}-{args.right}.json.gz"
)

timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
print(f"\nJob {job_name} started {timestamp}")
atoms_list: list[Atoms] = read(DataFiles.phonondb_pbe_103_structures.path, index=":")
atoms_list = atoms_list[args.left:args.right]  # Use the specified range

run_params = {
    "timestamp": timestamp,
    "model_name": model_name,
    "device": device,
    "versions": {dep: version(dep) for dep in ("numpy", "torch", "matbench_discovery")},
    "ase_optimizer": ase_optimizer,
    "cell_filter": "FrechetCellFilter",
    "max_steps": max_steps,
    "force_max": force_max,
    "symprec": symprec,
    "enforce_relax_symm": enforce_relax_symm,
    "conductivity_broken_symm": conductivity_broken_symm,
    "temperatures": temperatures,
    "displacement_distance": displacement_distance,
    "task_type": task_type,
    "job_name": job_name,
    "n_structures": len(atoms_list),
    "gpu": args.gpu,
    "left_index": args.left,
    "right_index": args.right,
}

with open(f"{out_dir}/run_params-gpu{args.gpu}-{args.left}-{args.right}.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Set up the relaxation and force set calculation
optim_cls: Callable[..., Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
force_results: dict[str, dict[str, Any]] = {}
kappa_results: dict[str, dict[str, Any]] = {}
tqdm_bar = tqdm(atoms_list, desc="Conductivity calculation: ", disable=not prog_bar)

for atoms in tqdm_bar:
    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    formula = atoms.get_chemical_formula()
    spg_num = MoyoDataset(MoyoAdapter.from_atoms(atoms)).number
    info_dict = {
        Key.desc: mat_id,
        Key.formula: formula,
        Key.spg_num: spg_num,
        "errors": [],
        "error_traceback": [],
    }

    tqdm_bar.set_postfix_str(mat_id, refresh=True)

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
                filtered_atoms, logfile=f"{out_dir}/relax_{mat_id}_gpu{args.gpu}.log"
            )
            optimizer.run(fmax=force_max, steps=max_steps)

            reached_max_steps = optimizer.nsteps >= max_steps
            if reached_max_steps:
                print(f"{mat_id=} reached {max_steps=} during relaxation")

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
        warnings.warn(f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2)
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
            fc2_supercell=atoms.info["fc2_supercell"],
            fc3_supercell=atoms.info["fc3_supercell"],
            q_point_mesh=atoms.info["q_point_mesh"],
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
                f"{mat_id=} has imaginary frequencies or broken symmetry", stacklevel=2
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
    force_out_path = f"{out_dir}/{today}-kappa-103-force-sets-gpu{args.gpu}-{args.left}-{args.right}.json.gz"
    df_force = pd.DataFrame(force_results).T
    df_force.index.name = Key.mat_id
    df_force.reset_index().to_json(force_out_path)
    print(f"Saved force results to {force_out_path}")
