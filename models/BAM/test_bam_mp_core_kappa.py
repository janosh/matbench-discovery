# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy>=1.23",
#   "pandas>=2.0",
#   "torch>=2.0",
#   "ase>=3.22",
#   "phono3py>=3.0",
#   "phonopy>=2.0",
#   "pymatviz>=0.15",
#   "tqdm>=4.65",
#   "matbench-discovery[phonons]",
#   "bam-torch",
# ]
# ///
"""Test BAM-MP-core on the matbench-discovery thermal-conductivity task.

This is the SevenNet kappa script with only the model/calculator setup adapted
for BAM-MP-core. Set ``BAM_MODEL_PKL`` to the BAM checkpoint path before running.
"""

from __future__ import annotations

import json
import os
import tempfile
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Literal

import pandas as pd
import torch
from ase.calculators.calculator import all_changes
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from bam_torch.tase.base_calculator import RACECalculator
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.enums import DataFiles
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from ase import Atoms
    from ase.optimize.optimize import Optimizer

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")


class BAMCalculator(RACECalculator):
    """RACECalculator wrapper exposing free_energy for ASE cell filters."""

    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "forces",
        "stress",
        "free_energy",
    ]

    def calculate(
        self,
        atoms: Atoms,
        properties: list[str] | tuple[str, ...] = ("energy",),
        system_changes: list[str] | None = None,
    ) -> None:  # type: ignore[override]
        if system_changes is None:
            system_changes = all_changes
        super().calculate(atoms, list(properties), system_changes)
        if "energy" in self.results:
            self.results["free_energy"] = self.results["energy"]


def prepare_bam_checkpoint(model_pkl: Path) -> Path:
    """Strip DDP ``module.`` prefixes from old BAM checkpoints when present."""
    ckpt = torch.load(model_pkl, map_location="cpu", weights_only=False)
    if not (
        isinstance(ckpt, dict)
        and "params" in ckpt
        and any(key.startswith("module.") for key in ckpt["params"])
    ):
        return model_pkl

    ckpt = dict(ckpt)
    ckpt["params"] = {
        key.replace("module.", "", 1): value for key, value in ckpt["params"].items()
    }
    with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp:
        torch.save(ckpt, tmp.name)
        tmp_path = Path(tmp.name)
    print(f"Stripped DDP module prefix into temporary checkpoint: {tmp_path}")
    return tmp_path


model_name = "BAM-MP-core"
model_variant = "bam-mp-core-1.0v"
device = "cuda" if torch.cuda.is_available() else "cpu"

model_pkl = Path(os.environ.get("BAM_MODEL_PKL", "model.pkl")).expanduser()
if not model_pkl.is_file():
    raise FileNotFoundError(
        "BAM checkpoint not found: "
        f"{model_pkl}. Set BAM_MODEL_PKL to the model.pkl path."
    )
model_pkl = prepare_bam_checkpoint(model_pkl.resolve())
calc = BAMCalculator(model=str(model_pkl), device=device)

# Relaxation parameters. These params follow SevenNet/EquiformerV3 kappa runs.
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
out_dir = os.environ.get("BAM_KAPPA_OUT_DIR", "./kappa_results")
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/{job_name}.json.gz"

timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
structures_path = Path(
    os.environ.get("BAM_KAPPA_STRUCTURES", DataFiles.phonondb_pbe_103_structures.path)
)
atoms_list = read(structures_path, index=":")

debug_n = int(os.environ.get("BAM_KAPPA_DEBUG", "0"))
if debug_n:
    atoms_list = atoms_list[:debug_n]

run_params = {
    "timestamp": timestamp,
    "model_name": model_name,
    "model_variant": model_variant,
    "model_pkl": str(model_pkl),
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

with open(f"{out_dir}/run_params.json", mode="w", encoding="utf-8") as file:
    json.dump(run_params, file, indent=4)

# Set up the relaxation and force set calculation
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
force_results: dict[str, dict[str, Any]] = {}
kappa_results: dict[str, dict[str, Any]] = {}
tqdm_bar = tqdm(atoms_list, desc="Conductivity calculation: ", disable=not prog_bar)

for atoms in tqdm_bar:
    mat_id = atoms.info[Key.mat_id]
    if "q_point_mesh" not in atoms.info and "q_mesh" in atoms.info:
        atoms.info["q_point_mesh"] = atoms.info["q_mesh"]
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

            log_path = f"{out_dir}/relax_{mat_id}.log"
            optimizer = optim_cls(filtered_atoms, logfile=log_path)  # ty: ignore[invalid-argument-type]
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

    except (ValueError, RuntimeError, OSError, KeyError) as exc:
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

    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        err_dict["errors"].append(f"ForceConstantError: {exc!r}")
        err_dict["error_traceback"].append(traceback.format_exc())
        kappa_results[mat_id] = info_dict | relax_dict | err_dict
        continue

    try:  # Calculate thermal conductivity
        ph3, kappa_dict, _ = ltc.calculate_conductivity(ph3, temperatures=temperatures)
        print(f"Calculated kappa for {mat_id}: {kappa_dict}")
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
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
