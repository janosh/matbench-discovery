"""Module with code for calculating the thermal conductivity of a material.

Code is adapted from https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/conductivity.py.
All credit to Balázs Póta, Paramvir Ahlawat, Gábor Csányi, Michele Simoncelli. See
https://arxiv.org/abs/2408.00755 for details.
It was ported to this repo in https://github.com/janosh/matbench-discovery/pull/196 to
implement parallelization across input structures which allows scaling thermal
conductivity metric to larger test sets.
"""

import warnings
from collections.abc import Sequence
from copy import deepcopy
from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.enums import MbdKey


def calculate_fc2_set(
    ph3: Phono3py, calculator: Calculator, pbar_kwargs: dict[str, Any] | None = None
) -> np.ndarray:
    """Calculate 2nd order force constants.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """
    print(f"Computing FC2 force set in {ph3.unitcell.formula}.")

    forces = []
    n_atoms = len(ph3.phonon_supercell)

    displacements = ph3.phonon_supercells_with_displacements
    for supercell in tqdm(
        displacements, desc=f"FC2 calculation: {ph3.unitcell.formula}", **pbar_kwargs
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
    **kwargs: Any,
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
        **kwargs (Any): Passed to Phono3py constructor.

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
        **kwargs,
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
    **kwargs: Any,
) -> tuple[Phono3py, dict[str, np.ndarray], Any]:
    """Calculate thermal conductivity.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate conductivity
        temperatures (list[float]): Temperatures to compute conductivity at in Kelvin.
        boundary_mfp (float): Mean free path in micrometer to calculate simple boundary
            scattering contribution to thermal conductivity. Defaults to 1e6.
        mode_kappa_thresh (float): Threshold for mode kappa consistency check. Defaults
            to 1e-6.
        **kwargs (Any): Passed to Phono3py.run_thermal_conductivity().

    Returns:
        tuple[Phono3py, dict[str, np.ndarray], Any]: (Phono3py object, conductivity
            dict, conductivity object)
    """
    ph3.init_phph_interaction(symmetrize_fc3q=False)

    ph3.run_thermal_conductivity(
        temperatures=temperatures,
        is_isotope=True,
        # use type="wigner" to include both wave-like coherence (kappa_c) and
        # particle-like (kappa_p) conductivity contributions
        conductivity_type="wigner",
        boundary_mfp=boundary_mfp,
        **kwargs,
    )

    kappa = ph3.thermal_conductivity

    kappa_dict = {
        MbdKey.kappa_tot_rta: deepcopy(kappa.kappa_TOT_RTA[0]),
        MbdKey.kappa_p_rta: deepcopy(kappa.kappa_P_RTA[0]),
        MbdKey.kappa_c: deepcopy(kappa.kappa_C[0]),
        Key.mode_weights: deepcopy(kappa.grid_weights),
        Key.q_points: deepcopy(kappa.qpoints),
        Key.ph_freqs: deepcopy(kappa.frequencies),
    }
    mode_kappa_total = kappa_dict[MbdKey.mode_kappa_tot_rta] = calc_mode_kappa_tot(
        deepcopy(kappa.mode_kappa_P_RTA[0]),
        deepcopy(kappa.mode_kappa_C[0]),
        deepcopy(kappa.mode_heat_capacities),
    )

    sum_mode_kappa_tot = mode_kappa_total.sum(
        axis=tuple(range(1, mode_kappa_total.ndim - 1))
    ) / np.sum(kappa_dict[Key.mode_weights])

    kappa_p_rta = kappa_dict[MbdKey.kappa_p_rta]
    if np.any((sum_mode_kappa_tot - kappa_p_rta) > mode_kappa_thresh):
        warnings.warn(
            f"Total mode kappa does not sum to total kappa. {sum_mode_kappa_tot=}, "
            f"{kappa_p_rta=}",
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
