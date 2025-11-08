"""Shared kappa calculation function for different ML potential models."""

import os
import traceback
import warnings
from collections.abc import Callable
from copy import deepcopy
from typing import TYPE_CHECKING, Any

import ase.optimize
import ase.optimize.sciopt
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, Filter, FrechetCellFilter
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatgen.core.structure import Structure
from pymatviz.enums import Key

from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from ase.optimize.optimize import Optimizer


def calc_kappa_for_structure(
    *,
    atoms: Atoms,
    calculator: Calculator,
    displacement_distance: float,
    temperatures: list[float],
    ase_optimizer: str,
    max_steps: int,
    force_max: float,
    symprec: float,
    enforce_relax_symm: bool,
    save_forces: bool,
    out_dir: str,
    task_id: int,
    ase_filter: str | None = None,
    conductivity_broken_symm: bool = False,
    ignore_imaginary_freqs: bool = False,
    formula_getter: Callable[[Atoms], str] | None = None,
) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    """Calculate thermal conductivity (kappa) for a single structure.

    This is a shared implementation used by different ML potential models
    (NequIP, Allegro, MACE, etc.).

    Args:
        atoms (Atoms): ASE Atoms object with fc2_supercell, fc3_supercell,
            q_point_mesh keys in its info dict.
        calculator (Calculator): ASE calculator to use for force calculations
        displacement_distance (float): Displacement distance for phono3py (Å)
        temperatures (list[float]): Temperatures in Kelvin for conductivity calculation
        ase_optimizer (str): ASE optimizer name (e.g., 'FIRE', 'BFGS', 'LBFGS')
        max_steps (int): Maximum relaxation steps
        force_max (float): Maximum force convergence criterion (eV/Å)
        symprec (float): Symmetry precision for spglib
        enforce_relax_symm (bool): Whether to enforce symmetry during relaxation
        save_forces (bool): Whether to save force sets to disk
        out_dir (str): Output directory for results
        task_id (int): Task ID for logging
        ase_filter (str | None): Cell filter for relaxation ('frechet' or 'exp').
            If None, uses FrechetCellFilter. Default None.
        conductivity_broken_symm (bool): Whether to calculate kappa if symmetry
            breaks. Only used if ignore_imaginary_freqs=False. Default False.
        ignore_imaginary_freqs (bool): Whether to ignore imaginary frequencies
            and calculate kappa anyway. Default False.
        formula_getter (Callable[[Atoms], str] | None): Custom function to extract
            formula from atoms. If None, uses atoms.get_chemical_formula().
            Default None.

    Returns:
        tuple[str, dict[str, Any], dict[str, Any] | None]:
            material ID, results dict, force results dict
    """
    formula = formula_getter(atoms) if formula_getter else atoms.get_chemical_formula()

    print(f"Calculating {Structure.from_ase_atoms(atoms).reduced_formula}")

    # Ensure arrays are writable
    atoms.arrays = {key: val.copy() for key, val in atoms.arrays.items()}

    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    info_dict: dict[str, Any] = {
        str(Key.mat_id): mat_id,
        str(Key.formula): formula,
    }
    err_dict: dict[str, list[str]] = {"errors": [], "error_traceback": []}

    # Select filter class
    if ase_filter in {"frechet", "exp"}:
        filter_cls: type[Filter] = {
            "frechet": FrechetCellFilter,
            "exp": ExpCellFilter,
        }[ase_filter]
    else:
        # Default to FrechetCellFilter if not specified (for MACE compatibility)
        filter_cls = FrechetCellFilter

    # Select optimizer class
    optimizer_dict = {
        "GPMin": ase.optimize.GPMin,
        "GOQN": ase.optimize.GoodOldQuasiNewton,
        "BFGSLineSearch": ase.optimize.BFGSLineSearch,
        "QuasiNewton": ase.optimize.BFGSLineSearch,
        "SciPyFminBFGS": ase.optimize.sciopt.SciPyFminBFGS,
        "BFGS": ase.optimize.BFGS,
        "LBFGSLineSearch": ase.optimize.LBFGSLineSearch,
        "SciPyFminCG": ase.optimize.sciopt.SciPyFminCG,
        "FIRE2": ase.optimize.FIRE2,
        "FIRE": ase.optimize.FIRE,
        "LBFGS": ase.optimize.LBFGS,
    }
    optim_cls: type[Optimizer] = optimizer_dict[ase_optimizer]

    # Initialize variables that might be needed in error handling
    relax_dict: dict[str, Any] = {
        "max_stress": None,
        "reached_max_steps": False,
        "broken_symmetry": False,
    }
    force_results = None
    # initial space group for symmetry breaking detection
    init_spg_num = MoyoDataset(MoyoAdapter.from_atoms(atoms), symprec=symprec).number

    # Relaxation
    try:
        atoms.calc = calculator
        if max_steps > 0:
            if enforce_relax_symm:
                atoms.set_constraint(FixSymmetry(atoms))
                filtered_atoms = filter_cls(atoms, mask=[True] * 3 + [False] * 3)
            else:
                filtered_atoms = filter_cls(atoms)

            os.makedirs(relax_dir := f"{out_dir}/relaxations", exist_ok=True)
            optimizer = optim_cls(filtered_atoms, logfile=f"{relax_dir}/{task_id}.log")
            optimizer.run(fmax=force_max, steps=max_steps)

            step_count = getattr(optimizer, "nsteps", None)  # Get optimizer step count
            if step_count is None:  # fallback to extract from state_dict if available
                state = getattr(optimizer, "state_dict", dict)()
                step_count = state.get("step", 0)

            reached_max_steps = step_count >= max_steps
            if reached_max_steps:
                print(f"Material {mat_id=} reached {max_steps=} during relaxation")

            # maximum residual stress component in for xx,yy,zz and xy,yz,xz
            # components separately result is a array of 2 elements
            max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

            atoms.calc = None
            atoms.constraints = None
            atoms.info = init_info | atoms.info

            # Check if symmetry was broken during relaxation
            moyo_cell = MoyoAdapter.from_atoms(atoms)
            relaxed_spg_num = MoyoDataset(moyo_cell, symprec=symprec).number
            broken_symmetry = init_spg_num != relaxed_spg_num

            relax_dict = {
                "max_stress": max_stress,
                "reached_max_steps": reached_max_steps,
                "broken_symmetry": broken_symmetry,
                "relaxed_space_group_number": relaxed_spg_num,
            }

    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        warnings.warn(f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        err_dict["errors"] += [f"RelaxError: {exc!r}"]
        err_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict | err_dict, None

    # Calculation of force sets
    try:
        ph3 = ltc.init_phono3py(
            atoms,
            fc2_supercell=atoms.info["fc2_supercell"],
            fc3_supercell=atoms.info["fc3_supercell"],
            q_point_mesh=atoms.info["q_point_mesh"],
            displacement_distance=displacement_distance,
            symprec=symprec,
        )

        ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
            ph3, calculator=calculator, pbar_kwargs={"disable": True}
        )

        # Lazy import to avoid circular dependency
        from matbench_discovery.phonons import check_imaginary_freqs

        has_imaginary_freqs = check_imaginary_freqs(freqs)
        freqs_dict = {Key.has_imag_ph_modes: has_imaginary_freqs, Key.ph_freqs: freqs}

        # Determine if conductivity calculation should proceed
        if ignore_imaginary_freqs:
            # NequIP/Allegro mode: ignore imaginary frequencies
            ltc_condition = True
        else:
            # MACE mode: check both imaginary freqs and broken symmetry
            broken_symmetry = relax_dict.get("broken_symmetry", False)
            ltc_condition = not has_imaginary_freqs and (
                not broken_symmetry or conductivity_broken_symm
            )

        if ltc_condition:
            fc3_set = ltc.calculate_fc3_set(
                ph3, calculator=calculator, pbar_kwargs={"position": task_id}
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            reason = []
            if has_imaginary_freqs:
                reason.append("imaginary frequencies")
            if relax_dict.get("broken_symmetry") and not conductivity_broken_symm:
                reason.append("broken symmetry")
            warnings.warn(
                f"{' and '.join(reason).capitalize()} detected for {mat_id}, "
                f"skipping FC3 and LTC calculation!",
                stacklevel=2,
            )
            fc3_set = []

        force_results = (
            {"fc2_set": fc2_set, "fc3_set": fc3_set} if save_forces else None
        )

        if not ltc_condition:
            return mat_id, info_dict | relax_dict | freqs_dict | err_dict, force_results

    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        err_dict["errors"] += [f"ForceConstantError: {exc!r}"]
        err_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict | err_dict, force_results

    # Calculation of conductivity
    try:
        ph3, kappa_dict, _cond = ltc.calculate_conductivity(
            ph3, temperatures=temperatures
        )
        return (
            mat_id,
            info_dict | relax_dict | freqs_dict | kappa_dict | err_dict,
            force_results,
        )

    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        err_dict["errors"] += [f"ConductivityError: {exc!r}"]
        err_dict["error_traceback"] += [traceback.format_exc()]
        return mat_id, info_dict | relax_dict | freqs_dict | err_dict, force_results
