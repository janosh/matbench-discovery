"""
Copyright (c) Meta Platforms, Inc. and its affiliates.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""
from __future__ import annotations
import glob
from pathlib import Path
from copy import deepcopy
import time
import json

import warnings
from typing import Any
import traceback

import random
import numpy as np
import torch
import pandas as pd
from tqdm import tqdm

from ase.io import read
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE

from k_srme import (
    aseatoms2str, 
    two_stage_relax, 
    glob2df, 
    ID, 
    STRUCTURES, 
    NO_TILT_MASK, 
    DFT_NONAC_REF
)
from k_srme.benchmark import get_metrics, process_benchmark_descriptors
from k_srme.utils import symm_name_map, get_spacegroup_number, check_imaginary_freqs
from k_srme.conductivity import (
    init_phono3py,
    get_fc2_and_freqs,
    get_fc3,
    calculate_conductivity,
)

from pymatviz.enums import Key
from matbench_discovery.enums import MbdKey

from fairchem.core import OCPCalculator

def seed_everywhere(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    
class kSRMERunner:
    
    def __init__(
        self, 
        seed,
        model_dir,
        save_name,
        identifier,
        atom_disp,
        num_jobs,
    ) -> None:
        self.seed = seed
        self.model_dir = model_dir
        self.save_name = save_name
        self.identifier = identifier
        self.atom_disp = atom_disp
        self.num_jobs = num_jobs
    
    def run(self, job_number=0) -> None:
        if_two_stage_relax = True  # Use two-stage relaxation enforcing symmetries
        max_steps = 300
        force_max = 1e-4  # Run until the forces are smaller than this in eV/A
        symprec = 1e-5
        enforce_relax_symm = True
        conductivity_broken_symm = False
        prog_bar = True
        save_forces = False  # Save force sets to
        filter_cls = FrechetCellFilter
        optim_cls = FIRE
        
        seed_everywhere(self.seed)
        
        model_ckpt = str(Path(self.model_dir) / "checkpoint.pt")
        calculator = OCPCalculator(checkpoint_path=model_ckpt, cpu=False, seed=0)
        calculator.trainer.scaler = None
        
        force_results: dict[str, dict[str, Any]] = {}
        kappa_results: dict[str, dict[str, Any]] = {}
        atoms_list = read(STRUCTURES, format="extxyz", index=":")
        
        if self.num_jobs > 0:
            atoms_list = atoms_list[job_number::self.num_jobs]
                
        tqdm_bar = tqdm(atoms_list, desc="Conductivity calculation: ", disable=not prog_bar)

        test_metrics = {}            
        save_dir = Path(self.model_dir) / self.save_name / f'{self.identifier}_{self.seed}'
        (save_dir).mkdir(parents=True, exist_ok=True)

        start_time = time.time()
        for atoms in tqdm_bar:
            mat_id = atoms.info[ID]
            init_info = deepcopy(atoms.info)
            mat_name = atoms.info["name"]
            mat_desc = f"{mat_name}-{symm_name_map[atoms.info['symm.no']]}"
            info_dict = {
                "desc": mat_desc,
                "name": mat_name,
                "initial_space_group_number": atoms.info["symm.no"],
                "errors": [],
                "error_traceback": [],
            }

            tqdm_bar.set_postfix_str(mat_desc, refresh=True)

            # Relaxation
            try:
                atoms.calc = calculator
                if max_steps > 0:
                    if not if_two_stage_relax:
                        if enforce_relax_symm:
                            atoms.set_constraint(FixSymmetry(atoms))
                            filtered_atoms = filter_cls(atoms, mask=NO_TILT_MASK)
                        else:
                            filtered_atoms = filter_cls(atoms)

                        optimizer = optim_cls(
                            filtered_atoms, logfile=save_dir / f"relax_{job_number}.log"
                        )
                        optimizer.run(fmax=force_max, steps=max_steps)

                        reached_max_steps = False
                        if optimizer.step == max_steps:
                            reached_max_steps = True
                            print(
                                f"Material {mat_desc=}, {mat_id=} reached max step {max_steps=} during relaxation."
                            )

                        # maximum residual stress component in for xx,yy,zz and xy,yz,xz components separately
                        # result is a array of 2 elements
                        max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

                        atoms.calc = None
                        atoms.constraints = None
                        atoms.info = init_info | atoms.info

                        symm_no = get_spacegroup_number(atoms, symprec=symprec)

                        relax_dict = {
                            "structure": aseatoms2str(atoms),
                            "max_stress": max_stress,
                            "reached_max_steps": reached_max_steps,
                            "relaxed_space_group_number": symm_no,
                            "broken_symmetry": symm_no
                            != init_info["initial_space_group_number"],
                        }

                    else:
                        atoms, relax_dict = two_stage_relax(
                            atoms,
                            fmax_stage1=force_max,
                            fmax_stage2=force_max,
                            steps_stage1=max_steps,
                            steps_stage2=max_steps,
                            Optimizer=optim_cls,
                            Filter=filter_cls,
                            allow_tilt=False,
                            log=save_dir / f"relax_{job_number}.log",
                            enforce_symmetry=enforce_relax_symm,
                        )

                        atoms.calc = None

            except Exception as exc:
                warnings.warn(f"Failed to relax {mat_name=}, {mat_id=}: {exc!r}")
                traceback.print_exc()
                info_dict["errors"].append(f"RelaxError: {exc!r}")
                info_dict["error_traceback"].append(traceback.format_exc())
                kappa_results[mat_id] = info_dict
                continue

            # Calculation of force sets
            try:
                ph3 = init_phono3py(atoms, log=False, symprec=symprec, displacement_distance=self.atom_disp)

                ph3, fc2_set, freqs = get_fc2_and_freqs(
                    ph3,
                    calculator=calculator,
                    log=False,
                    pbar_kwargs={"leave": False, "disable": not prog_bar},
                )

                imaginary_freqs = check_imaginary_freqs(freqs)
                freqs_dict = {"imaginary_freqs": imaginary_freqs, "frequencies": freqs}

                # if conductivity condition is met, calculate fc3
                ltc_condition = not imaginary_freqs and (
                    not relax_dict["broken_symmetry"] or conductivity_broken_symm #pylint: disable=E0606
                )

                if ltc_condition:
                    ph3, fc3_set = get_fc3(
                        ph3,
                        calculator=calculator,
                        log=False,
                        pbar_kwargs={"leave": False, "disable": not prog_bar},
                    )

                else:
                    fc3_set = []

                if save_forces:
                    force_results[mat_id] = {"fc2_set": fc2_set, "fc3_set": fc3_set}

                if not ltc_condition:
                    kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
                    warnings.warn(f"Material {mat_desc}, {mat_id} has imaginary frequencies.")
                    continue

            except Exception as exc:
                warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}")
                traceback.print_exc()
                info_dict["errors"].append(f"ForceConstantError: {exc!r}")
                info_dict["error_traceback"].append(traceback.format_exc())
                kappa_results[mat_id] = info_dict | relax_dict
                continue

            # Calculation of conductivity
            try:
                ph3, kappa_dict = calculate_conductivity(ph3, log=False)

            except Exception as exc:
                warnings.warn(f"Failed to calculate conductivity {mat_id}: {exc!r}")
                traceback.print_exc()
                info_dict["errors"].append(f"ConductivityError: {exc!r}")
                info_dict["error_traceback"].append(traceback.format_exc())
                kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
                continue

            kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict

        elapsed = time.time() - start_time
        test_metrics['running_time'] = elapsed
        
        df_kappa = pd.DataFrame(kappa_results).T
        df_kappa.index.name = ID
        df_kappa.reset_index().to_json(save_dir / f"conductivity_{self.num_jobs}-{job_number}.json.gz")
        
        # compute metrics
        pattern = str(save_dir / "conductivity_*-*.json.gz")
        file_list = list(glob.glob(pattern))
        if len(file_list) == self.num_jobs:
            df_mlp_results = glob2df(pattern, max_files=None).set_index(ID)
    
            # df_mlp_results = df_kappa.set_index(ID)
            df_dft_results = pd.read_json(DFT_NONAC_REF).set_index(ID)
            
            # assert len(df_mlp_results) == len(df_dft_results)
            
            df_mlp_filtered = df_mlp_results[df_mlp_results.index.isin(df_dft_results.index)]
            df_mlp_filtered = df_mlp_filtered.reindex(df_dft_results.index)
            df_mlp_processed = process_benchmark_descriptors(df_mlp_filtered, df_dft_results)
            mSRE, mSRME, rmseSRE, rmseSRME = get_metrics(df_mlp_filtered)

            test_metrics['mSRE'] = mSRE
            test_metrics['mSRME'] = mSRME
            test_metrics['rmseSRE'] = rmseSRE
            test_metrics['rmseSRME'] = rmseSRME

            df_mlp_processed.round(5)
            df_mlp_processed.index.name = ID
            df_mlp_processed.reset_index().to_json(str(save_dir / "k_srme.json.gz"))

            with open(save_dir / 'test_metrics.json', 'w', encoding='utf-8') as f:
                json.dump(test_metrics, f, ensure_ascii=False, indent=4)
                
            df_mlp_processed = df_mlp_processed.rename(columns={
                'mp_id': Key.mat_id,
                'sre': Key.sre,
                "srme": Key.srme,
                "weights": Key.mode_weights,
                "kappa_P_RTA": MbdKey.kappa_p_rta,
                "kappa_C": MbdKey.kappa_c,
                "kappa_TOT_RTA": MbdKey.kappa_tot_rta,
                "kappa_TOT_ave": MbdKey.kappa_tot_avg,
                "mode_kappa_TOT_ave": MbdKey.mode_kappa_tot_avg,
            })
            
            df_mlp_processed.to_json(str(save_dir / f"2025-03-17-kappa-103-FIRE-dist={self.atom_disp}-fmax=1e-4-symprec=1e-5.json.gz"))