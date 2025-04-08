import json
import os
import traceback
import warnings
from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from typing import Any, Literal

import pandas as pd
import torch
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from ase.spacegroup.symmetrize import check_symmetry
from orb_models.forcefield.calculator import ORBCalculator
from orb_models.forcefield.pretrained import ORB_PRETRAINED_MODELS
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from pymatviz.enums import Key
from tqdm import tqdm

import matbench_discovery.phonons.thermal_conductivity as ltc
from matbench_discovery import phonons
from matbench_discovery.enums import DataFiles
from matbench_discovery.metrics.phonons import calc_kappa_metrics_from_dfs

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")
warnings.filterwarnings("ignore", category=FutureWarning, module="torch")

# Model configuration
model_name = "orb-v3"
model_variant = "orb-v3-conservative-inf-mpa"  # ORB v3 model to evaluate
device = "cuda" if torch.cuda.is_available() else "cpu"
max_num_neighbors = 120

# Relaxation parameters
ase_optimizer: Literal["FIRE", "LBFGS"] = "FIRE"
max_steps = 300
force_max = 1e-4  # In eV/Å
symprec = 1e-5
displacement_distance = 0.03  # Displacement distance for phono3py
enforce_relax_symm = True
ignore_broken_symm = False
ignore_imaginary_freqs = False
is_plusminus = True
temperatures = [300]  # Temperatures to calculate conductivity at in Kelvin
save_forces = True  # Save force sets to file
deterministic = False
precision = "float64"

task_type = "LTC"  # lattice thermal conductivity
job_name = (
    f"{model_name}-phononDB-{task_type}-{ase_optimizer}_force{force_max}_sym{symprec}"
)
out_dir = f"./kappa_results_{model_variant}"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/{job_name}.json.gz"
force_sets_path = f"{out_dir}/force-sets.json.gz"

timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
atoms_list = read(DataFiles.phonondb_pbe_103_structures.path, index=":")

# Limit to only 10 structures
atoms_list = atoms_list[:10]

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
    "ignore_broken_symm": ignore_broken_symm,
    "ignore_imaginary_freqs": ignore_imaginary_freqs,
    "temperatures": temperatures,
    "displacement_distance": displacement_distance,
    "is_plusminus": is_plusminus,
    "task_type": task_type,
    "job_name": job_name,
    "n_structures": len(atoms_list),
    "max_num_neighbors": max_num_neighbors,
    "precision": precision,
    "deterministic": deterministic,
}

with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

print(f"Results will be saved to {out_dir}")
print(f"Using {device=}")

# Load and configure the ORB model
model = ORB_PRETRAINED_MODELS[model_variant]()
model.to(device)
calc = ORBCalculator(
    model,
    max_num_neighbors=max_num_neighbors,
    device=device,
)

if deterministic:
    torch.use_deterministic_algorithms(mode=True)


# Set up the optimizer class from string
optim_cls: Callable[..., Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

force_results: dict[str, dict[str, Any]] = {}
kappa_results: dict[str, dict[str, Any]] = {}
prog_bar = True  # Enable progress bar

tqdm_bar = tqdm(
    enumerate(atoms_list), desc="Conductivity calculation: ", disable=not prog_bar
)

for idx, atoms in tqdm_bar:
    # Use the same ID field as in original script
    mat_id = atoms.info.get(Key.mat_id, f"id-{len(kappa_results)}")
    init_info = deepcopy(atoms.info)
    formula = atoms.info.get("name", "unknown")

    tqdm_bar.set_postfix_str(mat_id, refresh=True)

    # Initialize info dictionary with material details
    info_dict = {
        "name": formula,
        "errors": [],
        "error_traceback": [],
    }

    # Initialize variables that might be needed in error handling
    relax_dict = {"max_stress": None, "reached_max_steps": False}
    force_results_item = None

    try:
        # Relaxation phase
        atoms.calc = calc
        if enforce_relax_symm:
            atoms.set_constraint(FixSymmetry(atoms))
            filtered_atoms = FrechetCellFilter(atoms, mask=[True] * 3 + [False] * 3)
        else:
            filtered_atoms = FrechetCellFilter(atoms)

        optimizer = optim_cls(
            filtered_atoms,
            logfile=f"{out_dir}/relax_{idx}.log",
        )

        pre_sym_group = check_symmetry(atoms, symprec).number
        optimizer.run(fmax=force_max, steps=max_steps)
        post_sym_group = check_symmetry(atoms, symprec).number

        reached_max_steps = optimizer.nsteps >= max_steps
        if reached_max_steps:
            print(f"Material {mat_id=} reached {max_steps=} during relaxation")

        # Maximum residual stress component
        max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

        atoms.calc = None
        atoms.constraints = None
        atoms.info = init_info | atoms.info

        relax_dict = {
            "max_stress": max_stress,
            "reached_max_steps": reached_max_steps,
            "broken_symmetry": pre_sym_group != post_sym_group,
        }

        if not ignore_broken_symm and pre_sym_group != post_sym_group:
            raise ValueError(
                f"Symmetry group changed from {pre_sym_group} to {post_sym_group}"
            )

    except Exception as exc:
        warnings.warn(f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"RelaxError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict
        continue

    # Force constants calculation
    try:
        unit_cell = PhonopyAtoms(
            atoms.symbols, cell=atoms.cell, positions=atoms.positions
        )
        ph3 = Phono3py(
            unitcell=unit_cell,
            supercell_matrix=atoms.info["fc3_supercell"],
            phonon_supercell_matrix=atoms.info["fc2_supercell"],
            primitive_matrix="auto",
            symprec=symprec,
        )
        ph3.mesh_numbers = atoms.info["q_mesh"]

        # Generate displacements in both positive and negative direction even if
        # symmetrically equivalent (different from other models!)
        ph3.generate_displacements(
            distance=displacement_distance, is_plusminus=is_plusminus
        )
        # Calculate force constants and frequencies
        ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
            ph3, calculator=calc, pbar_kwargs={"disable": True}
        )

        # Check for imaginary frequencies
        has_imaginary_freqs = phonons.check_imaginary_freqs(freqs)
        freqs_dict = {
            Key.has_imag_ph_modes: has_imaginary_freqs,
            Key.ph_freqs: freqs,
        }

        # Determine if we should continue calculating conductivity
        continue_computing_conductivity = (
            not has_imaginary_freqs or ignore_imaginary_freqs
        )

        if continue_computing_conductivity:
            fc3_set = ltc.calculate_fc3_set(
                ph3,
                calculator=calc,
                pbar_kwargs={"position": idx},
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            fc3_set = []

        if save_forces:
            force_results_item = {"fc2_set": fc2_set, "fc3_set": fc3_set}

        if not continue_computing_conductivity:
            kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
            warnings.warn(
                f"Skipping {mat_id} due to imaginary frequencies", stacklevel=2
            )
            if force_results_item is not None:
                force_results[mat_id] = force_results_item
            continue

    except Exception as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"ForceConstantError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict
        continue

    # Thermal conductivity calculation
    try:
        ph3, kappa_dict, _cond = ltc.calculate_conductivity(
            ph3, temperatures=temperatures
        )
        kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict
        if force_results_item is not None:
            force_results[mat_id] = force_results_item
    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        info_dict["errors"].append(f"ConductivityError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
        if force_results_item is not None:
            force_results[mat_id] = force_results_item

# Save results
df_kappa = pd.DataFrame(kappa_results).T
df_kappa.index.name = Key.mat_id
df_kappa.reset_index().to_json(out_path)
print(f"Saved kappa results to {out_path}")

if save_forces and force_results:
    df_force = pd.DataFrame(force_results).T
    df_force = pd.concat([df_kappa[[]].copy(), df_force], axis=1)
    df_force.index.name = Key.mat_id
    df_force.reset_index().to_json(force_sets_path)
    print(f"Saved force sets to {force_sets_path}")

try:
    print("Computing metrics against reference data...")
    df_dft = pd.read_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path).set_index(
        Key.mat_id
    )
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
