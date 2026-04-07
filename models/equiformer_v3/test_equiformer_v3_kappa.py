from __future__ import annotations

import json
import os
import random
import time
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import torch
import typer
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE
from fairchem.core import OCPCalculator
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.enums import DataFiles
from matbench_discovery.metrics import phonons
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc


def seed_everywhere(seed: int) -> None:
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


class KappaSRMERunner:
    def __init__(
        self,
        checkpoint_path,
        output_dir,
    ) -> None:

        self.checkpoint_path = str(checkpoint_path)
        self.output_dir = str(output_dir)

        # torch.set_default_dtype(torch.float64)

    def run(self) -> None:
        """
        Relaxation parameters
        1.  eSEN uses `max_steps` == 300 and `force_max` == 1e-4
        2.  Nequix uses `max_steps` == 500 and `force_max` == 0.02
        """
        max_steps = 300
        force_max = 1e-4  # Run until the forces are smaller than this in eV/A
        symprec = 1e-5
        enforce_relax_symm = True
        conductivity_broken_symm = False
        prog_bar = True
        save_forces = True  # Save force sets to file
        temperatures = [300]  # Temperatures to calculate conductivity at in Kelvin
        displacement_distance = 0.03  # Displacement distance for phono3py

        eval_every = 1

        seed_everywhere(0)

        # Setup model and calculator
        calculator = OCPCalculator(
            checkpoint_path=self.checkpoint_path, cpu=False, seed=0
        )
        calculator.trainer.scaler = None

        force_results: dict[str, dict[str, Any]] = {}
        kappa_results: dict[str, dict[str, Any]] = {}
        atoms_list = read(
            DataFiles.phonondb_pbe_103_structures.path, format="extxyz", index=":"
        )

        tqdm_bar = tqdm(
            atoms_list, desc="Conductivity calculation: ", disable=not prog_bar
        )

        test_metrics = {}
        save_dir = Path(self.output_dir)
        (save_dir).mkdir(parents=True, exist_ok=True)

        # Log run parameters
        timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
        run_params = {
            "timestamp": timestamp,
            "checkpoint_path": str(self.checkpoint_path),
            "device": "cuda" if torch.cuda.is_available() else "cpu",
            "versions": {
                dep: version(dep) for dep in ("numpy", "torch", "matbench_discovery")
            },
            "ase_optimizer": "FIRE",
            "cell_filter": "FrechetCellFilter",
            "max_steps": max_steps,
            "force_max": force_max,
            "symprec": symprec,
            "enforce_relax_symm": enforce_relax_symm,
            "conductivity_broken_symm": conductivity_broken_symm,
            "temperatures": temperatures,
            "displacement_distance": displacement_distance,
            "n_structures": len(atoms_list),
        }
        with open(save_dir / "run_params.json", mode="w") as file:
            json.dump(run_params, file, indent=4)

        start_time = time.time()
        index = 0
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

            # Relaxation - using standard approach from other scripts
            try:
                atoms.calc = calculator
                if max_steps > 0:
                    if enforce_relax_symm:
                        atoms.set_constraint(FixSymmetry(atoms))

                    # Use standard mask for no-tilt constraint
                    filtered_atoms = FrechetCellFilter(
                        atoms, mask=[True] * 3 + [False] * 3
                    )

                    optimizer = FIRE(
                        filtered_atoms, logfile=save_dir / f"relax_{mat_id}.log"
                    )
                    optimizer.run(fmax=force_max, steps=max_steps)

                    reached_max_steps = optimizer.nsteps >= max_steps
                    if reached_max_steps:
                        print(f"{mat_id=} reached {max_steps=} during relaxation")

                    max_stress = (
                        atoms.get_stress().reshape((2, 3), order="C").max(axis=1)
                    )
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
                warnings.warn(
                    f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2
                )
                traceback.print_exc()
                info_dict["errors"].append(f"RelaxError: {exc!r}")
                info_dict["error_traceback"].append(traceback.format_exc())
                kappa_results[mat_id] = info_dict | relax_dict
                with open(save_dir / "error.txt", "a") as txt_file:
                    txt_file.write(f"{mat_id}, relax\n")
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
                    calculator=calculator,
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
                        calculator=calculator,
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
                        f"{mat_id=} has imaginary frequencies or broken symmetry",
                        stacklevel=2,
                    )
                    continue

            except Exception as exc:
                warnings.warn(
                    f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2
                )
                traceback.print_exc()
                info_dict["errors"].append(f"ForceConstantError: {exc!r}")
                info_dict["error_traceback"].append(traceback.format_exc())
                kappa_results[mat_id] = info_dict | relax_dict
                with open(save_dir / "error.txt", "a") as txt_file:
                    txt_file.write(f"{mat_id}, forces\n")
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
                with open(save_dir / "error.txt", "a") as txt_file:
                    txt_file.write(f"{mat_id}, conductivity\n")
                continue

            kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict

            index = index + 1
            if index % eval_every == 0:
                save_kappa_results(kappa_results, save_dir)
                evaluate_tednet_kappa(f"{save_dir}/conductivity.json.gz")

        elapsed = time.time() - start_time
        test_metrics["running_time"] = elapsed

        # Save results
        save_kappa_results(kappa_results, save_dir)

        if save_forces:
            force_out_path = f"{save_dir}/force_sets.json.gz"
            df_force = pd.DataFrame(force_results).T
            df_force.index.name = Key.mat_id
            df_force.reset_index().to_json(force_out_path)
            print(f"Saved force results to {force_out_path}")


def save_kappa_results(kappa_results, save_dir):
    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    json_path = f"{save_dir}/conductivity.json.gz"
    df_kappa.reset_index().to_json(json_path)
    print(f"Saved kappa results to {json_path}")


def evaluate_tednet_kappa(path):
    pred_file = path
    df_ml = pd.read_json(pred_file).set_index("material_id")

    json_path = DataFiles.phonondb_pbe_103_kappa_no_nac.path
    df_dft = pd.read_json(json_path).set_index("material_id")

    df_ml_metrics = phonons.calc_kappa_metrics_from_dfs(df_ml, df_dft)

    kappa_srme_list = list(df_ml_metrics[Key.srme])
    kappa_sre = df_ml_metrics[Key.sre].mean()
    kappa_srme = df_ml_metrics[Key.srme].mean()

    dir = os.path.dirname(path)
    output_path = os.path.join(dir, "results.txt")
    with open(output_path, "a") as file:
        file.write(f"kappa_srme for each structure: {kappa_srme_list}\n")
        file.write(f"kappa_srme: \t{kappa_srme}\n")
        file.write(f"kappa_sre: \t{kappa_sre}\n")
    print(f"kappa_srme for each structure: {kappa_srme_list}")
    print(f"kappa_srme: \t{kappa_srme}")
    print(f"kappa_sre: \t{kappa_sre}")


def kappa_run_relaxation(
    checkpoint_path: Annotated[Path, typer.Option()],
    output_dir: Annotated[
        Path, typer.Option(help="Output path to write results files")
    ],
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    args = (
        checkpoint_path,
        output_dir,
    )
    runner = KappaSRMERunner(*args)
    runner.run()
    evaluate_tednet_kappa(f"{output_dir}/conductivity.json.gz")


if __name__ == "__main__":
    typer.run(kappa_run_relaxation)
