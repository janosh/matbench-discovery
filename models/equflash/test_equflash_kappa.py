# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
# "torch==2.8.0+cu126 ",
# "torch-geometric==2.6.1",
# "numpy==1.26.0",
# "scikit-learn==1.7.2",
# "spglib==2.6.0",
# "e3nn==0.5.6",
# "ase==3.26.0",
# "pymatgen==2025.10.7",
# "pymatviz>=0.16.0",
# "phono3py"==3.19.3,
# "flashTP_e3nn==0.1.0",
# "fairchem-core==1.10.0",
# "lmdb==1.6.2",
# "submitit==1.5.3",
# "matbench-discovery"==1.3.1,
# ]
#
# [tool.uv.sources]
# flashTP_e3nn = { git = "https://github.com/SNU-ARC/flashTP" }
# matbench-discovery = { path = "../../", editable = true }
# ///
import argparse
import datetime
import json
import os
import traceback
import warnings
from copy import deepcopy
from importlib.metadata import version
from typing import TYPE_CHECKING, Any, Literal, cast

import ase
import numpy as np
import pandas as pd
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, Filter, FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.spacegroup import get_spacegroup
from ase.utils import atoms_to_spglib_cell
from GGNN.common.calculator import UCalculator
from k_srme.conductivity import calculate_conductivity
from spglib import get_symmetry_dataset
from thermal_conductivity import get_fc3_batch
from tqdm import tqdm

from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons.thermal_conductivity import (
    get_fc2_and_freqs,
    init_phono3py,
)

if TYPE_CHECKING:
    from ase.optimize.optimize import Optimizer

ID = "mp_id"
NO_TILT_MASK = [True, True, True, False, False, False]
SYMM_NAME_MAP = {225: "rs", 186: "wz", 216: "zb"}


def log_symmetry(atoms: Atoms, symprec: float) -> Any:
    """Get symmetry dataset from atoms using spglib."""
    return get_symmetry_dataset(atoms_to_spglib_cell(atoms), symprec=symprec)


def two_stage_relax(
    atoms: Atoms,
    calculator: Calculator | None = None,
    fmax_stage1: float = 1e-4,
    fmax_stage2: float = 1e-4,
    steps_stage1: int = 300,
    steps_stage2: int = 300,
    *,
    enforce_symmetry: bool = True,
    symprec: float = 1e-5,
    allow_tilt: bool = False,
    optimizer: type["Optimizer"] = LBFGS,
    filter_ase: type[ase.filters.Filter] = FrechetCellFilter,
    symprec_tests: list[float] | None = None,
) -> tuple[Atoms, dict[str, Any]]:
    """Two-stage relaxation enforcing symmetry in first stage."""
    if calculator is not None:
        atoms.calc = calculator
    elif atoms.calc is None:
        raise ValueError("Atoms object does not have a calculator assigned")
    if symprec_tests is None:
        symprec_tests = [1e-5, 1e-4, 1e-3, 1e-1]

    _filter_kwargs = {}
    _optim_kwargs = {}

    if "name" in atoms.info:
        mat_name = atoms.info["name"]
    else:
        mat_name = f"{atoms.get_chemical_formula(mode='metal', empirical=True)}"
        f"-{get_spacegroup(atoms, symprec=symprec).no}"

    tilt_mask = None
    if not allow_tilt:
        tilt_mask = NO_TILT_MASK

    sym_init = log_symmetry(atoms, symprec)

    atoms.set_constraint(FixSymmetry(atoms))

    total_filter = filter_ase(atoms, mask=tilt_mask, **_filter_kwargs)
    dyn_stage1 = optimizer(total_filter, **_optim_kwargs)
    dyn_stage1.run(fmax=fmax_stage1, steps=steps_stage1)
    sym_stage1 = log_symmetry(atoms, symprec)

    if sym_stage1["number"] != sym_init["number"]:
        warnings.warn(
            f"Symmetry is not kept during FixSymmetry "
            f"relaxation of material {mat_name} in folder {os.getcwd()}",
            stacklevel=2,
        )

    max_stress_stage1 = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

    atoms_stage1 = atoms.copy()
    atoms.constraints = None

    dyn_stage2 = optimizer(total_filter, **_optim_kwargs)
    dyn_stage2.run(fmax=fmax_stage2, steps=steps_stage2)
    sym_stage2 = log_symmetry(atoms, symprec)

    sym_tests = {}
    if sym_init.number != sym_stage2.number:
        for symprec_test in symprec_tests:
            dataset_tests = log_symmetry(atoms, symprec_test)
            sym_tests[symprec_test] = dataset_tests.number

    max_stress_stage2 = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)
    if sym_stage1.number != sym_stage2.number and enforce_symmetry:
        redirected_to_symm = True
        atoms = atoms_stage1
        max_stress = max_stress_stage1
        sym_final = sym_stage1
        warnings.warn(
            f"Symmetry not kept after removing FixSymmetry constraint for "
            f"{mat_name} in {os.getcwd()}, redirecting to stage1",
            stacklevel=2,
        )
    else:
        redirected_to_symm = False
        sym_final = sym_stage2
        max_stress = max_stress_stage2

    reached_max_steps = (
        dyn_stage1.step == steps_stage1 or dyn_stage2.step == steps_stage2
    )

    relax_dict = {
        "max_stress": max_stress,
        "reached_max_steps": reached_max_steps,
        "relaxed_space_group_number": sym_final.number,
        "broken_symmetry": sym_final.number != sym_init.number,
        "symprec_tests": sym_tests,
        "redirected_to_symm": redirected_to_symm,
    }

    return_atoms = atoms

    return return_atoms, relax_dict


def get_spacegroup_number(atoms: Atoms, symprec: float = 1e-5) -> int:
    """Get space group number from atoms."""
    dataset = get_symmetry_dataset(atoms_to_spglib_cell(atoms), symprec=symprec)
    return dataset.number


def main() -> None:
    """Run thermal conductivity calculations with EquFlash model."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", "-c", required=True, help="Model checkpoint")
    parser.add_argument("--outdir", "-o", required=True, help="Output directory")
    parser.add_argument("--displacement", type=float, default=0.03)
    parser.add_argument("--rank", type=int, default=0)
    parser.add_argument("--worldsize", type=int, default=1)
    parser.add_argument("--structures", type=str, required=True)

    args = parser.parse_args()

    checkpoint = args.checkpoint
    device = "cuda"
    calc = UCalculator(checkpoint_path=checkpoint, cpu=False)

    ase_optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE"
    ase_filter: Literal["frechet", "exp"] = "frechet"
    if_two_stage_relax = True
    max_steps = 300
    force_max = 1e-4
    symprec = 1e-5
    enforce_relax_symm = True
    conductivity_broken_symm = False
    prog_bar = True
    save_forces = True
    task_type = "LTC"
    out_dir = args.outdir
    os.makedirs(out_dir, exist_ok=True)

    out_path = f"{out_dir}/conductivity-{args.rank:>03}.json.gz"

    timestamp = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S")
    struct_data_path = args.structures
    atoms_list = cast("list[Atoms]", read(struct_data_path, format="extxyz", index=":"))

    run_params = {
        "timestamp": timestamp,
        "k_srme_version": version("k_srme"),
        "checkpoint": checkpoint,
        "device": device,
        "versions": {dep: version(dep) for dep in ("numpy", "torch")},
        "ase_optimizer": ase_optimizer,
        "ase_filter": ase_filter,
        "if_two_stage_relax": if_two_stage_relax,
        "max_steps": max_steps,
        "force_max": force_max,
        "symprec": symprec,
        "enforce_relax_symm": enforce_relax_symm,
        "conductivity_broken_symm": conductivity_broken_symm,
        "slurm_array_task_count": args.worldsize,
        "task_type": task_type,
        "struct_data_path": os.path.basename(struct_data_path),
        "n_structures": len(atoms_list),
    }

    if args.rank == 0:
        with open(f"{out_dir}/run_params.json", "w") as f:
            json.dump(run_params, f, indent=4)

    atoms_list = atoms_list[args.rank :: args.worldsize]

    # Set up the relaxation and force set calculation
    filter_cls: type[Filter] = {
        "frechet": FrechetCellFilter,
        "exp": ExpCellFilter,
    }[ase_filter]
    optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

    force_results: dict[str, dict[str, Any]] = {}
    kappa_results: dict[str, dict[str, Any]] = {}

    print(f"{len(atoms_list)} on {args.rank}")
    for atoms in tqdm(atoms_list):
        mat_id = atoms.info[ID]
        init_info = deepcopy(atoms.info)
        mat_name = atoms.info["name"]
        mat_desc = f"{mat_name}-{SYMM_NAME_MAP[atoms.info['symm.no']]}"
        info_dict = {
            "desc": mat_desc,
            "name": mat_name,
            "initial_space_group_number": atoms.info["symm.no"],
            "errors": [],
            "error_traceback": [],
        }

        atoms.calc = calc
        if max_steps > 0:
            if not if_two_stage_relax:
                if enforce_relax_symm:
                    atoms.set_constraint(FixSymmetry(atoms))
                    filtered_atoms = filter_cls(atoms, mask=NO_TILT_MASK)
                else:
                    filtered_atoms = filter_cls(atoms)

                optimizer = optim_cls(
                    filtered_atoms, logfile=f"{out_dir}/relax_{args.rank}.log"
                )
                optimizer.run(fmax=force_max, steps=max_steps)

                reached_max_steps = False
                if optimizer.step == max_steps:
                    reached_max_steps = True
                    print(
                        f"Material {mat_desc=}, {mat_id=} reached max step "
                        f"{max_steps=} during relaxation."
                    )

                max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

                atoms.calc = None
                atoms.constraints = None
                atoms.info = init_info | atoms.info

                symm_no = get_spacegroup_number(atoms, symprec=symprec)

                relax_dict = {
                    "structure": atoms,
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
                    optimizer=optim_cls,
                    filter_ase=filter_cls,
                    allow_tilt=False,
                    enforce_symmetry=enforce_relax_symm,
                )
                relax_dict["ase_fc2_supercell"] = atoms.info["fc2_supercell"]
                relax_dict["ase_fc3_supercell"] = atoms.info["fc3_supercell"]
                relax_dict["ase_symbols"] = str(atoms.symbols)
                relax_dict["ase_cell"] = np.array(atoms.cell)
                relax_dict["ase_positions"] = atoms.positions
                relax_dict["ase_q_mesh"] = atoms.info["q_mesh"]

                atoms.calc = None

        info_dict["displacement"] = args.displacement
        try:
            ph3 = init_phono3py(
                atoms,
                fc2_supercell=atoms.info["fc2_supercell"],
                fc3_supercell=atoms.info["fc3_supercell"],
                q_point_mesh=atoms.info["q_mesh"],
                symprec=symprec,
                displacement_distance=args.displacement,
            )

            ph3, fc2_set, freqs = get_fc2_and_freqs(
                ph3,
                calculator=calc,
                pbar_kwargs={"leave": False, "disable": not prog_bar},
            )
            imaginary_freqs = check_imaginary_freqs(freqs, threshold=-1e-2)
            freqs_dict = {"imaginary_freqs": imaginary_freqs, "frequencies": freqs}

            ltc_condition = not imaginary_freqs and (
                not relax_dict["broken_symmetry"] or conductivity_broken_symm
            )

            if ltc_condition:
                ph3, fc3_set = get_fc3_batch(
                    ph3,
                    calculator=calc,
                    pbar_kwargs={"leave": False, "disable": not prog_bar},
                )
                ph3.forces = fc3_set
                ph3.produce_fc3(symmetrize_fc3r=True)
                ph3, kappa_dict = calculate_conductivity(ph3, log=False)
            else:
                fc3_set = []

            if save_forces:
                force_results[mat_id] = {"fc2_set": fc2_set, "fc3_set": fc3_set}

            if not ltc_condition:
                kappa_results[mat_id] = info_dict | relax_dict | freqs_dict
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

        kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | kappa_dict

    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = ID
    df_kappa.reset_index().to_json(out_path)

    df_force = pd.DataFrame(force_results).T
    df_force = pd.concat([df_kappa, df_force], axis=1)
    df_force.index.name = ID
    df_force.reset_index().to_json(out_path)


if __name__ == "__main__":
    main()
