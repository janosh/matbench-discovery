"""Module with code for calculating the thermal conductivity of a material.

The low-level phono3py force-constant and conductivity helpers are adapted from
https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/conductivity.py. All credit to
Balázs Póta, Paramvir Ahlawat, Gábor Csányi, Michele Simoncelli. See
https://arxiv.org/abs/2408.00755 for details. They were ported to this repo in
https://github.com/janosh/matbench-discovery/pull/196 to implement parallelization
across input structures which allows scaling thermal conductivity metric to larger test
sets. High-level orchestration lives in ``matbench_discovery.phonons.pipeline``.
"""

import warnings
from collections.abc import Callable, Sequence
from copy import deepcopy
from typing import TYPE_CHECKING, Literal

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms

from matbench_discovery.enums import Key

if TYPE_CHECKING:
    from phono3py.conductivity.calculators import LBTECalculator, RTACalculator


def ase_batch_forces(calculator: Calculator) -> Callable[[list[Atoms]], np.ndarray]:
    """Return a batch evaluator calling an ASE calculator on one structure at a time."""

    def evaluate_batch(batch_atoms: list[Atoms]) -> np.ndarray:
        for atoms in batch_atoms:
            atoms.calc = calculator
        return np.asarray([atoms.get_forces() for atoms in batch_atoms])

    return evaluate_batch


def batched_displacement_forces(
    displacements: Sequence[PhonopyAtoms | None],
    evaluate_batch: Callable[[list[Atoms]], np.ndarray],
    *,
    batch_size: int,
    n_atoms: int,
    max_evaluations: int | None = None,
) -> np.ndarray:
    """Evaluate the non-null displaced supercells, up to max_evaluations of them, in
    batches of batch_size periodic ASE structures, rejecting forces that are not finite
    real numbers of shape (batch length, n_atoms, 3).

    Returns:
        np.ndarray: (len(displacements), n_atoms, 3) forces, zero for null and
            unevaluated displacements.
    """
    if batch_size < 1 or n_atoms < 1:
        raise ValueError("batch_size and n_atoms must be positive")
    if max_evaluations is not None and max_evaluations < 0:
        raise ValueError("max_evaluations must be non-negative or None")
    entries = [
        (displacement_idx, supercell)
        for displacement_idx, supercell in enumerate(displacements)
        if supercell is not None
    ][:max_evaluations]
    forces = np.zeros((len(displacements), n_atoms, 3))
    for batch_start in range(0, len(entries), batch_size):
        batch = entries[batch_start : batch_start + batch_size]
        batch_atoms = [
            Atoms(cell.symbols, cell=cell.cell, positions=cell.positions, pbc=True)
            for _, cell in batch
        ]
        predicted = np.asarray(evaluate_batch(batch_atoms))
        expected_shape = (len(batch), n_atoms, 3)
        if (
            predicted.shape != expected_shape
            or not np.issubdtype(predicted.dtype, np.number)
            or np.iscomplexobj(predicted)
            or not np.isfinite(predicted).all()
        ):
            raise ValueError(
                f"Invalid displacement forces: shape={predicted.shape}, "
                f"dtype={predicted.dtype}, {expected_shape=}"
            )
        forces[[displacement_idx for displacement_idx, _ in batch]] = predicted
    return forces


def calculate_fc2_set(ph3: Phono3py, calculator: Calculator) -> np.ndarray:
    """Evaluate the FC2 displacement forces and assign them to ``ph3.phonon_forces``.

    Requires initializing Phono3py with an FC2 supercell matrix.
    """
    force_set = batched_displacement_forces(
        ph3.phonon_supercells_with_displacements,
        ase_batch_forces(calculator),
        batch_size=1,
        n_atoms=len(ph3.phonon_supercell),
    )
    ph3.phonon_forces = force_set
    return force_set


def calculate_fc3_set(
    ph3: Phono3py, calculator: Calculator, *, max_evaluations: int | None = None
) -> np.ndarray:
    """Evaluate the FC3 displacement forces and assign them to ``ph3.forces``.

    Only the first max_evaluations non-null displacements are evaluated (the rest
    stay zero to preserve the force-set shape); None evaluates all.
    """
    force_set = batched_displacement_forces(
        ph3.supercells_with_displacements,
        ase_batch_forces(calculator),
        batch_size=1,
        n_atoms=len(ph3.supercell),
        max_evaluations=max_evaluations,
    )
    ph3.forces = force_set
    return force_set


def init_phono3py(
    atoms: Atoms,
    *,
    fc2_supercell: np.ndarray,
    fc3_supercell: np.ndarray,
    q_point_mesh: tuple[int, int, int] = (20, 20, 20),
    displacement_distance: float = 0.03,
    symprec: float = 1e-5,
    is_plusminus: bool | Literal["auto"] = "auto",
) -> Phono3py:
    """Initialize Phono3py from ASE Atoms and generate its FC2/FC3 displacements.

    The default displacement_distance (Å) matches the PhononDB-PBE reference force
    sets, is_plusminus="auto" defers to phono3py. q_point_mesh must be 3 positive
    integers: whole floats are accepted, anything else raises instead of truncating.
    """
    mesh = np.asarray(q_point_mesh, dtype=float)
    if mesh.shape != (3,) or not all(np.isfinite(mesh) & (mesh >= 1) & (mesh % 1 == 0)):
        raise ValueError(f"q_point_mesh needs 3 positive integers, got {q_point_mesh}")
    unit_cell = PhonopyAtoms(atoms.symbols, cell=atoms.cell, positions=atoms.positions)
    ph3 = Phono3py(
        unitcell=unit_cell,
        supercell_matrix=fc3_supercell,
        phonon_supercell_matrix=fc2_supercell,
        primitive_matrix="auto",
        symprec=symprec,
    )
    ph3.mesh_numbers = mesh.astype(int).tolist()
    ph3.generate_displacements(
        distance=displacement_distance, is_plusminus=is_plusminus
    )
    return ph3


def fc2_frequencies(ph3: Phono3py, fc2_set: np.ndarray) -> np.ndarray:
    """Fit symmetrized FC2 to fc2_set and return the (n_bz_grid, n_bands) phonon
    frequencies on phono3py's BZ grid.
    """
    if ph3.mesh_numbers is None:
        raise ValueError("phono3py has no mesh_numbers to calculate phonons on")
    ph3.phonon_forces = fc2_set
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.init_phph_interaction(symmetrize_fc3q=False)
    ph3.run_phonon_solver()
    return ph3.get_phonon_data()[0]


def get_fc2_and_freqs(
    ph3: Phono3py, calculator: Calculator
) -> tuple[Phono3py, np.ndarray, np.ndarray]:
    """Return ph3, its FC2 force set from an ASE calculator and fc2_frequencies."""
    fc2_set = calculate_fc2_set(ph3, calculator)
    return ph3, fc2_set, fc2_frequencies(ph3, fc2_set)


def calculate_conductivity(
    ph3: Phono3py,
    temperatures: Sequence[float],
    boundary_mfp: float = 1e6,
    mode_kappa_thresh: float = 1e-6,
) -> tuple[Phono3py, dict[str, np.ndarray], RTACalculator | LBTECalculator]:
    """Calculate the thermal conductivity at temperatures (K) with isotope and
    boundary_mfp (µm) boundary scattering, warning if the total mode kappa deviates
    from the total kappa by more than mode_kappa_thresh.

    Returns:
        tuple: (ph3, conductivity dict, phono3py conductivity calculator)
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
        Key.kappa_tot_rta: deepcopy(extra["kappa_TOT_RTA"][0]),
        Key.kappa_p_rta: deepcopy(extra["kappa_P_RTA"][0]),
        Key.kappa_c: deepcopy(extra["kappa_C"][0]),
        Key.mode_weights: deepcopy(kappa.grid_weights),
        Key.q_points: deepcopy(kappa.qpoints),
        Key.ph_freqs: deepcopy(kappa.frequencies),
    }
    mode_kappa_total = kappa_dict[Key.mode_kappa_tot_rta] = calc_mode_kappa_tot(
        extra["mode_kappa_P_RTA"][0], extra["mode_kappa_C"][0], mode_cv
    )

    sum_mode_kappa_tot = mode_kappa_total.sum(
        axis=tuple(range(1, mode_kappa_total.ndim - 1))
    ) / np.sum(kappa_dict[Key.mode_weights])

    kappa_tot_rta = kappa_dict[Key.kappa_tot_rta]
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
    """Total (particle-like + wave-like) thermal conductivity per phonon mode.

    Args:
        mode_kappa_p_rta (np.ndarray): (T, q-points, bands, xyz) particle-like RTA
        mode_kappa_coherence (np.ndarray): (T, q-points, bands, bands, xyz) coherence
        heat_capacity (np.ndarray): (T, q-points, bands) mode heat capacities

    Returns:
        np.ndarray: (T, q-points, bands, xyz) total mode kappa
    """
    # Temporarily silence divide warnings since we handle NaN values below
    with np.errstate(divide="ignore", invalid="ignore"):
        mode_kappa_c_per_mode = 2 * (  # None equiv to np.newaxis
            (mode_kappa_coherence * heat_capacity[:, :, :, None, None])
            / (heat_capacity[:, :, :, None, None] + heat_capacity[:, :, None, :, None])
        ).sum(axis=2)
    mode_kappa_c_per_mode[np.isnan(mode_kappa_c_per_mode)] = 0
    return mode_kappa_c_per_mode + mode_kappa_p_rta
