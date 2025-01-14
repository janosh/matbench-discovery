"""This script runs SevenNet kappa predictions in parallel using Ray.
Can be run both locally and on a cluster with automatic resource scaling.
"""

import json
import os
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from typing import TYPE_CHECKING, Any, Literal

import pandas as pd
import ray
import torch
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from k_srme import ID, NO_TILT_MASK, STRUCTURES, aseatoms2str, two_stage_relax
from k_srme.conductivity import (
    calculate_conductivity,
    get_fc2_and_freqs,
    get_fc3,
    init_phono3py,
)
from k_srme.utils import check_imaginary_freqs, get_spacegroup_number, symm_name_map
from sevenn.sevennet_calculator import SevenNetCalculator
from tqdm import tqdm

from matbench_discovery import today

if TYPE_CHECKING:
    from ase.optimize.optimize import Optimizer

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")


@ray.remote
def process_structure(
    *,  # force keyword-only arguments
    atoms: Atoms,
    calc: SevenNetCalculator,
    ase_optimizer: str,
    ase_filter: str,
    if_two_stage_relax: bool,
    max_steps: int,
    force_max: float,
    symprec: float,
    enforce_relax_symm: bool,
    conductivity_broken_symm: bool,
    save_forces: bool,
    out_dir: str,
    task_id: int,
) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    """Process a single structure with SevenNet.

    Args:
        atoms: Input structure
        calc: SevenNet calculator
        ase_optimizer: ASE optimizer to use
        ase_filter: ASE filter to use
        if_two_stage_relax: Whether to use two-stage relaxation
        max_steps: Maximum number of optimization steps
        force_max: Maximum force tolerance
        symprec: Symmetry precision
        enforce_relax_symm: Whether to enforce symmetry during relaxation
        conductivity_broken_symm: Whether to calculate conductivity if symmetry broken
        save_forces: Whether to save force sets
        out_dir: Output directory
        task_id: Task ID for logging

    Returns:
        Tuple of (material ID, results dict, force results dict)
    """
    # Create a deep copy of the atoms object to avoid ray read-only issues
    atoms = atoms.copy()
    # Ensure arrays are writable
    atoms.arrays = {k: v.copy() for k, v in atoms.arrays.items()}

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

    filter_cls: type[ExpCellFilter | FrechetCellFilter] = {
        "frechet": FrechetCellFilter,
        "exp": ExpCellFilter,
    }[ase_filter]

    optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

    # Initialize variables that might be needed in error handling
    relax_dict = {
        "structure": aseatoms2str(atoms),
        "max_stress": None,
        "reached_max_steps": False,
        "relaxed_space_group_number": atoms.info["symm.no"],
        "broken_symmetry": False,
    }
    force_results = None

    # Relaxation
    try:
        atoms.calc = calc
        if max_steps > 0:
            if not if_two_stage_relax:
                if enforce_relax_symm:
                    atoms.set_constraint(FixSymmetry(atoms))
                    filtered_atoms = filter_cls(atoms, mask=NO_TILT_MASK)
                else:
                    filtered_atoms = filter_cls(atoms)

                optimizer = optim_cls(
                    filtered_atoms, logfile=f"{out_dir}/relax_{task_id}.log"
                )
                optimizer.run(fmax=force_max, steps=max_steps)

                reached_max_steps = optimizer.step == max_steps
                if reached_max_steps:
                    print(
                        f"Material {mat_desc=}, {mat_id=} reached {max_steps=} during "
                        "relaxation"
                    )

                # maximum residual stress component in for xx,yy,zz and xy,yz,xz
                # components separately result is a array of 2 elements
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
                    log=f"{out_dir}/relax_{task_id}.log",
                    enforce_symmetry=enforce_relax_symm,
                )

                atoms.calc = None

    except Exception as exc:
        warnings.warn(f"Failed to relax {mat_name=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"RelaxError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict, None

    # Calculation of force sets
    try:
        ph3 = init_phono3py(atoms, log=False, symprec=symprec)

        ph3, fc2_set, freqs = get_fc2_and_freqs(
            ph3,
            calculator=calc,
            log=False,
            pbar_kwargs={"leave": False, "disable": True},
        )

        imaginary_freqs = check_imaginary_freqs(freqs)
        freqs_dict = {"imaginary_freqs": imaginary_freqs, "frequencies": freqs}

        # if conductivity condition is met, calculate fc3
        ltc_condition = not imaginary_freqs and (
            not relax_dict["broken_symmetry"] or conductivity_broken_symm
        )

        if ltc_condition:
            ph3, fc3_set = get_fc3(
                ph3,
                calculator=calc,
                log=False,
                pbar_kwargs={"leave": False, "disable": True},
            )
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
        info_dict["errors"].append(f"ForceConstantError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict, force_results

    # Calculation of conductivity
    try:
        ph3, kappa_dict, _cond = calculate_conductivity(ph3, log=False)
        return mat_id, info_dict | relax_dict | freqs_dict | kappa_dict, force_results

    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        info_dict["errors"].append(f"ConductivityError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict | freqs_dict, force_results


def main() -> None:
    """Main function to run SevenNet kappa predictions."""
    # Model parameters
    module_dir = os.path.dirname(__file__)
    model_name = "SevenNet_l3i5"
    checkpoint = f"{module_dir}/sevennet-l3i5/checkpoint_l3i5.pth"
    if not os.path.isfile(checkpoint):
        raise FileNotFoundError("SevenNet checkpoint not found")

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Relaxation parameters
    ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
    ase_filter: Literal["frechet", "exp"] = "frechet"
    if_two_stage_relax = True  # Use two-stage relaxation enforcing symmetries
    max_steps = 300
    force_max = 1e-4  # Run until the forces are smaller than this in eV/A

    # Symmetry parameters
    symprec = (
        1e-5  # symmetry precision for enforcing relaxation and conductivity calculation
    )
    enforce_relax_symm = True  # Enforce symmetry with during relaxation if broken
    # Conductivity to be calculated if symmetry group changed during relaxation
    conductivity_broken_symm = False
    save_forces = True  # Save force sets to file

    # Ray initialization
    if os.getenv("RAY_HEAD_ADDRESS"):
        # Connect to existing Ray cluster
        ray.init(address=os.getenv("RAY_HEAD_ADDRESS"))
    else:
        # Start Ray locally with optimized settings for M3 Max
        ray.init(
            num_cpus=8,  # Use 8/14 cores (leaving some for system + efficiency cores)
            num_gpus=1,  # M3 Max GPU will be treated as 1 GPU
        )

    print(f"\nConnected to Ray cluster: {ray.cluster_resources()}")
    ray_mem = ray.available_resources().get("memory", 0) / 1e9
    print(f"Available memory: {ray_mem:.1f} GB")
    obj_store_mem = ray.available_resources().get("object_store_memory", 0)
    print(f"Object store memory: {obj_store_mem / 1e9:.1f} GB")

    task_type = "LTC"  # lattice thermal conductivity
    job_name = (
        f"{model_name}-phononDB-{task_type}-{ase_optimizer}"
        f"{'_2SR' if if_two_stage_relax else ''}_force{force_max}_sym{symprec}"
    )

    out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
    os.makedirs(out_dir, exist_ok=True)

    timestamp = f"{datetime.now().astimezone():%Y-%m-%d@%H-%M-%S}"
    struct_data_path = STRUCTURES
    print(f"\nJob {job_name} started {timestamp}")

    print(f"Read data from {struct_data_path}")
    atoms_list: list[Atoms] = read(struct_data_path, format="extxyz", index=":")

    # Save run parameters
    run_params = dict(
        model_name=model_name,
        checkpoint=checkpoint,
        device=device,
        ase_optimizer=ase_optimizer,
        ase_filter=ase_filter,
        if_two_stage_relax=if_two_stage_relax,
        max_steps=max_steps,
        force_max=force_max,
        symprec=symprec,
        enforce_relax_symm=enforce_relax_symm,
        conductivity_broken_symm=conductivity_broken_symm,
        task_type=task_type,
        struct_data_path=struct_data_path,
        n_structures=len(atoms_list),
        versions={dep: version(dep) for dep in ("numpy", "torch", "ray")},
    )

    with open(f"{out_dir}/run_params.json", mode="w") as file:
        json.dump(run_params, file, indent=4)

    # Create SevenNet calculator
    calc = SevenNetCalculator(model=checkpoint, device=device)

    # Process structures in parallel
    futures = [
        process_structure.remote(
            atoms=atoms,
            calc=calc,
            ase_optimizer=ase_optimizer,
            ase_filter=ase_filter,
            if_two_stage_relax=if_two_stage_relax,
            max_steps=max_steps,
            force_max=force_max,
            symprec=symprec,
            enforce_relax_symm=enforce_relax_symm,
            conductivity_broken_symm=conductivity_broken_symm,
            save_forces=save_forces,
            out_dir=out_dir,
            task_id=idx,
        )
        for idx, atoms in enumerate(atoms_list)
    ]

    # Process results as they complete
    kappa_results: dict[str, dict[str, Any]] = {}
    force_results: dict[str, dict[str, Any]] = {}

    print("\nProcessing structures...")
    for future in tqdm(futures, desc=f"Predicting kappa with {model_name}"):
        mat_id, result_dict, force_dict = ray.get(future)
        kappa_results[mat_id] = result_dict
        if force_dict is not None:
            force_results[mat_id] = force_dict

        # Save intermediate results
        df_kappa = pd.DataFrame(kappa_results).T
        df_kappa.index.name = ID
        df_kappa.reset_index().to_json(f"{out_dir}/conductivity.json.gz")

        if save_forces:
            df_force = pd.DataFrame(force_results).T
            df_force = pd.concat([df_kappa, df_force], axis=1)
            df_force.index.name = ID
            df_force.reset_index().to_json(f"{out_dir}/force_sets.json.gz")

    print(f"\nResults saved to {out_dir!r}")


if __name__ == "__main__":
    main()
