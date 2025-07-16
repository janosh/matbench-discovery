"""
Copyright (c) Meta Platforms, Inc. and its affiliates.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

from __future__ import annotations
import os
import glob
import json
import random
import time
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, Annotated
import typer
from matbench_discovery import timestamp, today
from submitit import AutoExecutor


import pandas as pd
import torch
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
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")


def seed_everywhere(seed: int) -> None:
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


class KappaSRMERunner:
    def __init__(
        self,
        seed: int,
        model_dir: str,
        out_dir: str,
        save_name: str,
        identifier: str,
        atom_disp: float,
        num_jobs: int,
    ) -> None:
        self.seed = seed
        self.model_dir = model_dir
        self.out_dir = out_dir
        self.save_name = save_name
        self.identifier = identifier
        self.atom_disp = atom_disp
        self.num_jobs = num_jobs

    def run(self, 
        env_vars: None, 
        job_number: int = 0
        ) -> None:

        # 设置 CUDA_VISIBLE_DEVICES 环境变量
        if env_vars:
            os.environ.update(env_vars)
            
        # Relaxation parameters
        max_steps = 300
        force_max = 1e-4  # Run until the forces are smaller than this in eV/A
        symprec = 1e-5
        enforce_relax_symm = True
        conductivity_broken_symm = False
        prog_bar = True
        save_forces = False  # Save force sets to file
        temperatures = [300]  # Temperatures to calculate conductivity at in Kelvin

        seed_everywhere(self.seed)

        # Setup model and calculator
        calculator = OCPCalculator(checkpoint_path=self.model_dir, cpu=False, seed=0)
        calculator.trainer.scaler = None

        force_results: dict[str, dict[str, Any]] = {}
        kappa_results: dict[str, dict[str, Any]] = {}
        atoms_list = read(
            DataFiles.phonondb_pbe_103_structures.path, format="extxyz", index=":"
        )

        if self.num_jobs > 0:
            atoms_list = atoms_list[job_number :: self.num_jobs]
        else:
            atoms_list = atoms_list[:3] # for debug
            
        tqdm_bar = tqdm(
            atoms_list, desc="Conductivity calculation: ", disable=not prog_bar
        )

        test_metrics = {}
        save_dir = self.out_dir

        # Log run parameters
        timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
        run_params = {
            "timestamp": timestamp,
            "model_dir": str(self.model_dir),
            "device": "cuda" if torch.cuda.is_available() else "cpu",
            "versions": {
                dep: version(dep) for dep in ("numpy", "torch", "matbench_discovery")
            },
            "optimizer": "FIRE",
            "max_steps": max_steps,
            "force_max": force_max,
            "symprec": symprec,
            "enforce_relax_symm": enforce_relax_symm,
            "conductivity_broken_symm": conductivity_broken_symm,
            "temperatures": temperatures,
            "displacement_distance": self.atom_disp,
            "n_structures": len(atoms_list),
        }
        with open(save_dir / "run_params.json", mode="w") as file:
            json.dump(run_params, file, indent=4)

        start_time = time.time()
        for atoms in tqdm_bar:
            mat_id = atoms.info.get(Key.mat_id, f"id-{len(kappa_results)}")
            init_info = deepcopy(atoms.info)
            formula = atoms.info.get("name", "unknown")
            spg_num = MoyoDataset(MoyoAdapter.from_atoms(atoms)).number
            info_dict = {
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

                    # reached_max_steps = optimizer.step >= max_steps
                    reached_max_steps = optimizer.nsteps >= max_steps
                    if reached_max_steps:
                        print(f"{mat_id=} reached {max_steps=} during relaxation.")

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
                continue

            # Calculation of force sets
            try:
                # Initialize phono3py with the relaxed structure
                ph3 = ltc.init_phono3py(
                    atoms,
                    fc2_supercell=atoms.info.get("fc2_supercell", [2, 2, 2]),
                    fc3_supercell=atoms.info.get("fc3_supercell", [2, 2, 2]),
                    q_point_mesh=atoms.info.get("q_point_mesh", [10, 10, 10]),
                    displacement_distance=self.atom_disp,
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

        elapsed = time.time() - start_time
        test_metrics["running_time"] = elapsed

        # Save results
        df_kappa = pd.DataFrame(kappa_results).T
        df_kappa.index.name = Key.mat_id
        json_path = f"{save_dir}/conductivity_{self.num_jobs}-{job_number}.json.gz"
        df_kappa.reset_index().to_json(json_path)
        print(f"Saved kappa results to {json_path}")

        if save_forces:
            force_out_path = (
                save_dir / f"force_sets_{self.num_jobs}-{job_number}.json.gz"
            )
            df_force = pd.DataFrame(force_results).T
            df_force.index.name = Key.mat_id
            df_force.reset_index().to_json(force_out_path)
            print(f"Saved force results to {force_out_path}")

        # Compute metrics if we've collected all results
        pattern = str(save_dir / "conductivity_*-*.json.gz")
        file_list = list(glob.glob(pattern))
        if len(file_list) == self.num_jobs:
            # Load all results
            all_dfs = []
            for file_path in file_list:
                try:
                    all_dfs.append(pd.read_json(file_path))
                except Exception as exc:
                    warnings.warn(f"Failed to read {file_path}: {exc!r}", stacklevel=2)

            if all_dfs:
                combined_df = pd.concat(all_dfs).reset_index()
                combined_df.to_json(save_dir / "kappa-103-FIRE.json.gz")
                print(
                    f"Combined results saved to {save_dir / 'kappa-103-FIRE.json.gz'}"
                )

                with open(
                    save_dir / "test_metrics.json", mode="w", encoding="utf-8"
                ) as file:
                    json.dump(test_metrics, file, ensure_ascii=False, indent=4)


def run_kappa(
    checkpoint_path: Annotated[Path, typer.Option(help="路径，包含 checkpoint.pt")],
    out_path: Annotated[Path, typer.Option(help="Output path to write results files")],
    model_name: Annotated[str, typer.Option(help="输出目录的子目录名")],
    identifier: Annotated[str, typer.Option(help="运行标识，用于文件夹命名")],
    atom_disp: Annotated[float, typer.Option(help="用于计算声子性质的位移距离")] = 0.03,
    seed: Annotated[int, typer.Option(help="随机种子")] = 42,
    num_jobs: Annotated[int, typer.Option(help="并行任务数")] = 10,
    debug: Annotated[bool, typer.Option(help="调试模式，仅运行一个 job")] = False,
    slurm_timeout: Annotated[int, typer.Option(help="slurm 超时时间（小时）")] = 80,
) -> None:
    # 设置输出路径
    base_dir = Path(out_path) / model_name / f"{identifier}_{seed}"
    base_dir.mkdir(parents=True, exist_ok=True)

    out_path = base_dir / timestamp
    os.makedirs(out_path)

    # 配置 submitit Executor
    executor = AutoExecutor(folder=out_path)
    executor.update_parameters(
        cpus_per_task=1,
        gpus_per_task=1,
        nodes=1,
        timeout_min=60 * slurm_timeout,
        mem_gb=64,
    )

    #继续填充函数
    if debug:
        # 调试模式，只运行一个任务，使用GPU 0
        env_vars = {"CUDA_VISIBLE_DEVICES": "0"}
        job = executor.submit(
            KappaSRMERunner(
                seed=seed,
                model_dir=checkpoint_path,
                out_dir=out_path,
                save_name=model_name,
                identifier=identifier,
                atom_disp=atom_disp,
                num_jobs=1,  # 设置为0处理前3个样本
            ).run,
            env_vars=env_vars,
            job_number=0,
        )
    else:
        # 正常模式，提交多个任务
        jobs = []
        with executor.batch():
            for job_number in range(num_jobs):
                # 为每个任务分配不同的 GPU
                env_vars = {"CUDA_VISIBLE_DEVICES": str(job_number)}
                job = executor.submit(
                    KappaSRMERunner(
                        seed=seed,
                        model_dir=checkpoint_path,
                        out_dir=out_path,
                        save_name=model_name,
                        identifier=identifier,
                        atom_disp=atom_disp,
                        num_jobs=num_jobs,
                    ).run,
                    env_vars=env_vars,
                    job_number=job_number,
                )
                jobs.append(job)
    
if __name__ == "__main__":
    typer.run(run_kappa)