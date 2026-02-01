import json
import os
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import torch
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from sevenn.calculator import SevenNetCalculator
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.enums import DataFiles
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")

model_name = "sevennet"
model_variant = "sevennet-omni-i12"
device = "cuda" if torch.cuda.is_available() else "cpu"

calc_kwargs: dict[str, Any] = {
    "sevennet-0": {"model": "7net-0"},
    "sevennet-l3i5": {"model": "7net-l3i5"},
    "sevennet-mf-ompa": {"model": "7net-mf-ompa", "modal": "mpa"},
    "sevennet-omni-i12": {"model": "7net-omni-i12", "modal": "mpa"},
}[model_variant]
calc_kwargs["device"] = device

# Will be removed after integrating model checkpoint download into sevenn package
checkpoint_urls = {
    "sevennet-omni-i12": "https://figshare.com/ndownloader/files/60977863",
}
if model_variant in checkpoint_urls:
    cache_dir = Path.home() / ".cache" / "sevennet"
    cache_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = cache_dir / f"checkpoint_{model_variant.replace('-', '_')}.pth"

    if not checkpoint_path.exists():
        print(f"Downloading {model_variant} checkpoint to {checkpoint_path}...")
        import requests

        response = requests.get(checkpoint_urls[model_variant], stream=True, timeout=30)
        response.raise_for_status()
        with open(checkpoint_path, "wb") as f:
            f.writelines(response.iter_content(chunk_size=8192))
        print("Download complete.")
    else:
        print(f"Using cached checkpoint: {checkpoint_path}")
    calc_kwargs["model"] = str(checkpoint_path)

calc = SevenNetCalculator(**calc_kwargs)

# Relaxation parameters. These params are for reproducing 7net-mf-ompa.
ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
max_steps = 300
force_max = 1e-4
symprec = 1e-5
enforce_relax_symm = True
conductivity_broken_symm = False
prog_bar = True
save_forces = True  # Save force sets to file
temperatures = [300]  # Temperatures to calculate conductivity at in Kelvin
displacement_distance = 0.03  # Displacement distance for phono3py

task_type = "LTC"  # lattice thermal conductivity
job_name = (
    f"{today}-kappa-103-{ase_optimizer}-dist={displacement_distance}-"
    f"fmax={force_max}-{symprec=}"
)
out_dir = "./kappa_results"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/{job_name}.json.gz"

timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
atoms_list = read(DataFiles.phonondb_pbe_103_structures.path, index=":")

run_params = {
    "timestamp": timestamp,
    "model_name": model_name,
    "model_variant": model_variant,
    "device": device,
    "versions": {dep: version(dep) for dep in ("numpy", "torch")},
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
}

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

# Set up the relaxation and force set calculation
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
force_results: dict[str, dict[str, Any]] = {}
kappa_results: dict[str, dict[str, Any]] = {}
tqdm_bar = tqdm(atoms_list, desc="Conductivity calculation: ", disable=not prog_bar)

for atoms in tqdm_bar:
    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    formula = atoms.get_chemical_formula()

    spg_num = MoyoDataset(MoyoAdapter.from_atoms(atoms)).number

    info_dict: dict[str, Any] = {
        str(Key.mat_id): mat_id,
        str(Key.formula): formula,
        str(Key.spg_num): spg_num,
    }
    err_dict: dict[str, list[str]] = {"errors": [], "error_traceback": []}

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
                filtered_atoms, logfile=f"{out_dir}/relax_{mat_id}.log"
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
        err_dict["errors"].append(f"RelaxError: {exc!r}")
        err_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | err_dict
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
            ph3, calculator=calc, pbar_kwargs={"leave": False, "disable": not prog_bar}
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
        err_dict["errors"].append(f"ForceConstantError: {exc!r}")
        err_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | err_dict
        continue

    try:  # Calculate thermal conductivity
        ph3, kappa_dict, _ = ltc.calculate_conductivity(ph3, temperatures=temperatures)
        print(f"Calculated kappa for {mat_id}: {kappa_dict}")
    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        err_dict["errors"].append(f"ConductivityError: {exc!r}")
        err_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | err_dict
        continue

    kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict | err_dict

# Save results
df_kappa = pd.DataFrame(kappa_results).T
df_kappa.index.name = Key.mat_id
df_kappa.reset_index().to_json(out_path)

if save_forces:
    force_out_path = f"{out_dir}/{today}-kappa-103-force-sets.json.gz"
    df_force = pd.DataFrame(force_results).T
    df_force.index.name = Key.mat_id
    df_force.reset_index().to_json(force_out_path)
