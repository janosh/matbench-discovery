"""Force-constant and lattice thermal-conductivity helpers for the kappa benchmark.

The protocol (finite displacements for fc2 in one supercell and fc3 in another, Wigner
transport in the relaxation-time approximation with isotope and boundary scattering)
follows https://github.com/MPA2suite/k_SRME (Póta, Ahlawat, Csányi, Simoncelli,
https://arxiv.org/abs/2408.00755), ported in
https://github.com/janosh/matbench-discovery/pull/196. All numerics run in matscikit;
high-level orchestration lives in ``matbench_discovery.phonons.pipeline``.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from matscikit import AnharmonicPhonons, Structure
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey

PlusMinus = bool | Literal["auto"]


@dataclass
class KappaCalculation:
    """Phonons of one material plus the q-point mesh its conductivity uses."""

    phonons: AnharmonicPhonons
    mesh: tuple[int, int, int]


def init_phonons(
    structure: Structure,
    *,
    fc2_supercell: np.ndarray,
    fc3_supercell: np.ndarray,
    q_point_mesh: tuple[int, int, int] = (20, 20, 20),
    displacement_distance: float = 0.03,
    symprec: float = 1e-5,
    is_plusminus: PlusMinus = "auto",
) -> KappaCalculation:
    """Set up fc2/fc3 supercells and their symmetry-reduced displacements.

    Args:
        structure (Structure): Unit cell.
        fc2_supercell (np.ndarray): Supercell matrix for 2nd order force constants.
        fc3_supercell (np.ndarray): Supercell matrix for 3rd order force constants.
        q_point_mesh (tuple[int, int, int]): Mesh size for q-point sampling. Defaults
            to (20, 20, 20).
        displacement_distance (float): Displacement distance in Å. Defaults to 0.03,
            matching the PhononDB-PBE reference force sets.
        symprec (float): Symmetry precision. Defaults to 1e-5.
        is_plusminus (bool | Literal["auto"]): Whether to add minus displacements
            (always, never, or when symmetry cannot supply them). Defaults to "auto".

    Returns:
        KappaCalculation: Phonons with generated displacements and the mesh.
    """

    def as_int_matrix(matrix: np.ndarray) -> np.ndarray:
        """Integer supercell matrix (metadata may store whole floats)."""
        values = np.asarray(matrix, dtype=float)
        if np.abs(values - np.rint(values)).max() > 1e-8:
            raise ValueError(f"supercell matrix must be integer, got {values.tolist()}")
        return np.rint(values).astype(int)

    phonons = AnharmonicPhonons(
        structure,
        as_int_matrix(fc3_supercell),
        phonon_supercell_matrix=as_int_matrix(fc2_supercell),
        primitive_matrix="auto",
        symprec=symprec,
    )
    plus_minus = None if is_plusminus == "auto" else bool(is_plusminus)
    phonons.generate_displacements(
        distance=displacement_distance, plus_minus=plus_minus
    )
    mesh_values = np.asarray(q_point_mesh, dtype=float).ravel()
    if (
        len(mesh_values) != 3
        or not np.isfinite(mesh_values).all()
        or (mesh_values < 1).any()
        or (mesh_values != np.rint(mesh_values)).any()
    ):
        raise ValueError(f"q_point_mesh needs 3 positive integers, got {q_point_mesh}")
    mesh = tuple(int(n) for n in mesh_values)
    return KappaCalculation(phonons, mesh)  # ty: ignore[invalid-argument-type]


def calculate_conductivity(
    calc: KappaCalculation,
    temperatures: Sequence[float],
    boundary_mfp: float = 1e6,
) -> dict[str, np.ndarray]:
    """Wigner-RTA lattice thermal conductivity (populations + coherences).

    Args:
        calc (KappaCalculation): Phonons with fc2 and fc3 set.
        temperatures (list[float]): Temperatures in Kelvin.
        boundary_mfp (float): Boundary mean free path in micrometer. Defaults to 1e6.

    Returns:
        dict[str, np.ndarray]: kappa tensors (temperature, 6) in W/(m K) and
            mode-resolved data on the irreducible q-points: mode_kappa_tot_rta
            (temperature, q-point, band, 6) adds to each mode's populations term its
            heat-capacity share of the band-pair coherences (k_SRME's split, summing
            to kappa_tot_rta times the grid size).
    """
    result = calc.phonons.thermal_conductivity(
        list(calc.mesh),
        [float(temperature) for temperature in temperatures],
        isotope=True,
        boundary_mfp_um=boundary_mfp,
    )
    return {
        MbdKey.kappa_tot_rta: result.kappa_tot_rta,
        MbdKey.kappa_p_rta: result.kappa_p_rta,
        MbdKey.kappa_c: result.kappa_c,
        MbdKey.mode_kappa_tot_rta: result.mode_kappa_tot_rta,
        Key.mode_weights: result.weights,
        Key.q_points: result.q_points,
        Key.ph_freqs: result.frequencies,
    }
