"""Module with code for calculating the thermal conductivity of a material.

The low-level phono3py force-constant and conductivity helpers are adapted from
https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/conductivity.py. All credit to
Balázs Póta, Paramvir Ahlawat, Gábor Csányi, Michele Simoncelli. See
https://arxiv.org/abs/2408.00755 for details. They were ported to this repo in
https://github.com/janosh/matbench-discovery/pull/196 to implement parallelization
across input structures which allows scaling thermal conductivity metric to larger test
sets. The high-level per-structure orchestrator ``calc_kappa_for_structure`` (relax ->
FC2/freqs -> FC3 -> conductivity) is matbench-discovery's own and shared across all
ML-potential kappa model scripts.
"""

import os
import traceback
import warnings
from collections.abc import Callable, Sequence
from copy import deepcopy
from typing import TYPE_CHECKING, Any

import ase.optimize
import ase.optimize.sciopt
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, Filter, FrechetCellFilter
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core.structure import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.enums import MbdKey

if TYPE_CHECKING:
    from ase.optimize.optimize import Optimizer
    from phono3py.conductivity.calculators import LBTECalculator, RTACalculator


def calculate_fc2_set(
    ph3: Phono3py, calculator: Calculator, pbar_kwargs: dict[str, Any] | None = None
) -> np.ndarray:
    """Calculate 2nd order force constants. Requires initializing Phono3py with an FC2
    supercell matrix.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """
    print(f"Computing FC2 force set in {ph3.unitcell.formula}.")

    forces: list[np.ndarray] = []
    n_atoms = len(ph3.phonon_supercell)

    displacements = ph3.phonon_supercells_with_displacements
    for supercell in tqdm(
        displacements,
        desc=f"FC2 calculation: {ph3.unitcell.formula}",
        **pbar_kwargs or {},
    ):
        if supercell is not None:
            atoms = Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            atoms.calc = calculator
            force = atoms.get_forces()
        else:
            force = np.zeros((n_atoms, 3))
        forces += [force]

    force_set = np.array(forces)
    ph3.phonon_forces = force_set
    return force_set


def calculate_fc3_set(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
) -> np.ndarray:
    """Calculate 3rd order force constants.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """
    forces: list[np.ndarray] = []
    n_atoms = len(ph3.supercell)

    desc = f"FC3 calculation: {ph3.unitcell.formula}"
    task_idx = (pbar_kwargs or {}).get("position")
    if task_idx:
        desc = f"{task_idx}. {desc}"
    displacements = ph3.supercells_with_displacements
    for supercell in tqdm(displacements, desc=desc, **pbar_kwargs or {}):
        if supercell is None:
            forces += [np.zeros((n_atoms, 3))]
        else:
            atoms = Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            atoms.calc = calculator
            forces += [atoms.get_forces()]

    force_set = np.array(forces)
    ph3.forces = force_set
    return force_set


def init_phono3py(
    atoms: Atoms,
    *,
    fc2_supercell: np.ndarray,
    fc3_supercell: np.ndarray,
    q_point_mesh: tuple[int, int, int] = (20, 20, 20),
    displacement_distance: float = 0.01,
    symprec: float = 1e-5,
) -> Phono3py:
    """Initialize Phono3py object from ASE Atoms.

    Args:
        atoms (Atoms): ASE Atoms object to initialize from.
        fc2_supercell (np.ndarray): Supercell matrix for 2nd order force constants.
        fc3_supercell (np.ndarray): Supercell matrix for 3rd order force constants.
        q_point_mesh (tuple[int, int, int]): Mesh size for q-point sampling. Defaults
            to (20, 20, 20).
        displacement_distance (float): Displacement distance for force calculations.
            Defaults to 0.01.
        symprec (float): Symmetry precision for finding space group. Defaults to 1e-5.

    Returns:
        Phono3py: Initialized Phono3py object

    Raises:
        ValueError: If required metadata is missing from atoms.info
    """
    unit_cell = PhonopyAtoms(atoms.symbols, cell=atoms.cell, positions=atoms.positions)
    ph3 = Phono3py(
        unitcell=unit_cell,
        supercell_matrix=fc3_supercell,
        phonon_supercell_matrix=fc2_supercell,
        primitive_matrix="auto",
        symprec=symprec,
    )
    ph3.mesh_numbers = q_point_mesh

    ph3.generate_displacements(distance=displacement_distance)

    return ph3


def get_fc2_and_freqs(
    ph3: Phono3py, calculator: Calculator, pbar_kwargs: dict[str, Any] | None = None
) -> tuple[Phono3py, np.ndarray, np.ndarray]:
    """Calculate 2nd order force constants and phonon frequencies.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        tuple[Phono3py, np.ndarray, np.ndarray]: Tuple of (Phono3py object, force
            constants array, frequencies array [shape: (n_bz_grid, n_bands)])

    Raises:
        ValueError: If mesh_numbers not set
    """
    if ph3.mesh_numbers is None:
        raise ValueError(
            "mesh_numbers was not found in phono3py object and was not provided as "
            "an argument when calculating phonons from phono3py object."
        )

    pbar_kwargs = {"leave": False} | (pbar_kwargs or {})
    fc2_set = calculate_fc2_set(ph3, calculator, pbar_kwargs=pbar_kwargs)

    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.init_phph_interaction(symmetrize_fc3q=False)
    ph3.run_phonon_solver()

    freqs, _eigvecs, _grid = ph3.get_phonon_data()

    return ph3, fc2_set, freqs


def load_force_sets(
    ph3: Phono3py, fc2_set: np.ndarray, fc3_set: np.ndarray
) -> Phono3py:
    """Load pre-computed force sets into Phono3py object.

    Args:
        ph3 (Phono3py): Phono3py object to load force sets into
        fc2_set (np.ndarray): 2nd order force constants array
        fc3_set (np.ndarray): 3rd order force constants array

    Returns:
        Phono3py: Phono3py object with loaded force sets
    """
    ph3.phonon_forces = fc2_set
    ph3.forces = fc3_set
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)

    return ph3


def calculate_conductivity(
    ph3: Phono3py,
    temperatures: Sequence[float],
    boundary_mfp: float = 1e6,
    mode_kappa_thresh: float = 1e-6,
) -> tuple[Phono3py, dict[str, np.ndarray], "RTACalculator | LBTECalculator"]:
    """Calculate thermal conductivity.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate conductivity
        temperatures (list[float]): Temperatures to compute conductivity at in Kelvin.
        boundary_mfp (float): Mean free path in micrometer to calculate simple boundary
            scattering contribution to thermal conductivity. Defaults to 1e6.
        mode_kappa_thresh (float): Threshold for mode kappa consistency check. Defaults
            to 1e-6.

    Returns:
        tuple[Phono3py, dict[str, np.ndarray], RTACalculator]: (Phono3py object,
            conductivity dict, conductivity object)
    """
    ph3.init_phph_interaction(symmetrize_fc3q=False)

    ph3.run_thermal_conductivity(
        temperatures=temperatures,
        is_isotope=True,
        # use MS-SMM19 (Wigner transport equation) to include both wave-like
        # coherence (kappa_c) and particle-like (kappa_p) conductivity contributions
        transport_type="MS-SMM19",
        boundary_mfp=boundary_mfp,
    )

    if (kappa := ph3.thermal_conductivity) is None:
        raise RuntimeError(f"thermal conductivity calculation failed for {ph3=}")
    extra, mode_cv = kappa.get_extra_kappa_output(), kappa.mode_heat_capacities
    if extra is None or mode_cv is None:
        raise RuntimeError(f"missing kappa output for {ph3=}")

    kappa_dict: dict[str, np.ndarray] = {
        MbdKey.kappa_tot_rta: deepcopy(extra["kappa_TOT_RTA"][0]),
        MbdKey.kappa_p_rta: deepcopy(extra["kappa_P_RTA"][0]),
        MbdKey.kappa_c: deepcopy(extra["kappa_C"][0]),
        Key.mode_weights: deepcopy(kappa.grid_weights),
        Key.q_points: deepcopy(kappa.qpoints),
        Key.ph_freqs: deepcopy(kappa.frequencies),
    }
    mode_kappa_total = kappa_dict[MbdKey.mode_kappa_tot_rta] = calc_mode_kappa_tot(
        deepcopy(extra["mode_kappa_P_RTA"][0]),
        deepcopy(extra["mode_kappa_C"][0]),
        deepcopy(mode_cv),
    )

    sum_mode_kappa_tot = mode_kappa_total.sum(
        axis=tuple(range(1, mode_kappa_total.ndim - 1))
    ) / np.sum(kappa_dict[Key.mode_weights])

    kappa_tot_rta = kappa_dict[MbdKey.kappa_tot_rta]
    if np.any(np.abs(sum_mode_kappa_tot - kappa_tot_rta) > mode_kappa_thresh):
        warnings.warn(
            f"Total mode kappa does not sum to total kappa. {sum_mode_kappa_tot=}, "
            f"{kappa_tot_rta=}",
            stacklevel=2,
        )

    return ph3, kappa_dict, kappa


def calc_mode_kappa_tot(
    mode_kappa_p_rta: np.ndarray,
    mode_kappa_coherence: np.ndarray,
    heat_capacity: np.ndarray,
) -> np.ndarray:
    """Calculate total mode kappa from particle-like RTA and coherence terms.

    Args:
        mode_kappa_p_rta (np.ndarray): Mode kappa from particle-like RTA with shape
            (T, q-points, bands, xyz)
        mode_kappa_coherence (np.ndarray): Mode kappa from wave-like coherence with
            shape (T, q-points, bands, bands, xyz)
        heat_capacity (np.ndarray): Mode heat capacities with shape
            (T, q-points, bands)

    Returns:
        np.ndarray: Total (particle-like + wave-like) thermal conductivity per phonon
            mode with shape (T, q-points, bands)
    """
    # Temporarily silence divide warnings since we handle NaN values below
    with np.errstate(divide="ignore", invalid="ignore"):
        mode_kappa_c_per_mode = 2 * (  # None equiv to np.newaxis
            (mode_kappa_coherence * heat_capacity[:, :, :, None, None])
            / (heat_capacity[:, :, :, None, None] + heat_capacity[:, :, None, :, None])
        ).sum(axis=2)

    mode_kappa_c_per_mode[np.isnan(mode_kappa_c_per_mode)] = 0

    return mode_kappa_c_per_mode + mode_kappa_p_rta


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
    **_kwargs: object,
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
        **_kwargs: Additional keywords (unused).

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
            optimizer = optim_cls(filtered_atoms, logfile=f"{relax_dir}/{task_id}.log")  # ty: ignore[invalid-argument-type]
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
        ph3 = init_phono3py(
            atoms,
            fc2_supercell=atoms.info["fc2_supercell"],
            fc3_supercell=atoms.info["fc3_supercell"],
            q_point_mesh=atoms.info["q_point_mesh"],
            displacement_distance=displacement_distance,
            symprec=symprec,
        )

        ph3, fc2_set, freqs = get_fc2_and_freqs(
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
            fc3_set = calculate_fc3_set(
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
        ph3, kappa_dict, _cond = calculate_conductivity(ph3, temperatures=temperatures)
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
