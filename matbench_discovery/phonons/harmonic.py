"""Harmonic phonon data derived from second-order force constants.

The kappa pipeline only needs mesh frequencies to gate on imaginary modes and feed
phono3py's conductivity solver, but the same FC2 also determines the full harmonic
picture of a material: band structure along the Brillouin-zone high-symmetry path with
eigenvectors and group velocities, the tetrahedron-method density of states on the
benchmark q-mesh, and the resulting harmonic thermal properties. Everything here is
cheap relative to the force evaluations, so the pipeline records it for every material
as a sidecar to the conductivity record. The site turns the eigenvectors into
animated mode displacements (matterviz ``PhononModeExplorer``), so their layout
follows phonopy's ``band.yaml`` convention: mass-weighted dynamical-matrix
eigenvectors normalized to ``sum |e|^2 = 1``, stored as ``[re, im]`` pairs per atom
and Cartesian direction. Eigenvectors are kept at every path q-point (tens of MB gzipped
per model for the 103 PhononDB primitives, on the order of the prediction file) rather
than only at the labeled points the site animates, so bands can later be tracked through
crossings by eigenvector overlap when comparing ML and DFT band structures.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    from phono3py.api_phono3py import Phono3py
    from phonopy import Phonopy

HARMONIC_SCHEMA_VERSION = 1
# seekpath's default q-point spacing along the path in 1/A (2 pi convention): the FCC Cu
# GAMMA-X segment gets 69 points, in the range of phonopy's 51-point band.conf default
BAND_PATH_REFERENCE_DISTANCE = 0.025
FREQUENCY_UNIT = "THz"
THERMAL_T_MIN, THERMAL_T_MAX, THERMAL_T_STEP = 0.0, 1000.0, 10.0
# eigenvector components are rounded to this many decimals: matterviz re-normalizes
# within 1e-3 of |e|^2 = 1 and the rounding keeps the sidecar a few MB per material
EIGENVECTOR_DECIMALS = 5
# seekpath's standardized primitive must be the same lattice as phonopy's primitive
# (possibly re-oriented and re-based) for the path's q-points to transfer
LATTICE_ATOL = 1e-4


def seekpath_band_path(phonon: Phonopy) -> dict[str, Any]:
    """Sample the seekpath high-symmetry path of a phonopy primitive cell.

    seekpath standardizes the cell (its own orientation and primitive basis, which
    can differ from phonopy's even for undistorted PhononDB cells, e.g. by a 4-fold
    rotation for zincblende), so the path's fractional coordinates are mapped back
    through Cartesian space into phonopy's reciprocal primitive basis exactly like
    seekpath's own ``get_path_orig_cell``.

    Returns:
        seekpath's explicit-path result with ``explicit_kpoints_rel`` rewritten to
        fractional q-points in phonopy's reciprocal primitive basis, plus
        ``explicit_kpoints_linearcoord``, ``explicit_kpoints_labels`` and
        ``explicit_segments`` (half-open index pairs; connected segments share their
        endpoint).
    """
    import seekpath

    path = seekpath.get_explicit_k_path(
        phonon.primitive.totuple(), reference_distance=BAND_PATH_REFERENCE_DISTANCE
    )
    own_lattice = np.asarray(phonon.primitive.cell)
    # Cartesian row vectors map standardized -> input orientation via `@ rotation`
    rotation = np.asarray(path["rotation_matrix"])
    # both cells describe the same lattice iff phonopy's basis vectors, rotated into the
    # standardized orientation, are unimodular integer combinations of seekpath's
    basis_change = own_lattice @ rotation.T @ np.linalg.inv(path["primitive_lattice"])
    integer_diff = np.abs(basis_change - np.round(basis_change)).max()
    if integer_diff > LATTICE_ATOL or round(abs(np.linalg.det(basis_change))) != 1:
        raise ValueError(
            "seekpath standardized primitive is not a re-based phonopy primitive "
            f"(basis change {basis_change.round(4).tolist()}); the band path "
            "q-coordinates would be mislabeled"
        )
    k_cartesian = (  # 1/A with 2 pi, in the input cell's orientation
        np.asarray(path["explicit_kpoints_rel"])
        @ path["reciprocal_primitive_lattice"]
        @ rotation
    )
    path["explicit_kpoints_rel"] = k_cartesian @ own_lattice.T / (2 * np.pi)
    return path


def phonopy_from_phono3py(ph3: Phono3py) -> Phonopy:
    """Build a Phonopy object sharing phono3py's cells, symmetry, and FC2."""
    from phonopy import Phonopy

    if ph3.fc2 is None:
        raise ValueError("phono3py has no fc2; run produce_fc2 before harmonic data")
    supercell_matrix = (
        ph3.phonon_supercell_matrix
        if ph3.phonon_supercell_matrix is not None
        else ph3.supercell_matrix
    )
    phonon = Phonopy(
        ph3.unitcell,
        supercell_matrix=supercell_matrix,
        primitive_matrix=ph3.primitive_matrix,
        symprec=ph3.symmetry.tolerance,
    )
    phonon.force_constants = ph3.fc2
    return phonon


def harmonic_phonon_data(phonon: Phonopy, *, mesh: Sequence[int]) -> dict[str, Any]:
    """Compute band structure, eigenvectors, DOS, and thermal properties from FC2.

    Args:
        phonon: Phonopy object with force constants set.
        mesh: q-point mesh for the tetrahedron DOS and thermal properties.

    Returns:
        JSON-ready mapping. Eigenvectors along the seekpath band path have shape
        ``(n_q, n_bands, n_atoms, 3, 2)`` in phonopy's ``[re, im]`` convention,
        rounded to ``EIGENVECTOR_DECIMALS`` after the normalization check.
    """
    path = seekpath_band_path(phonon)
    labels = path["explicit_kpoints_labels"]
    phonon.run_qpoints(
        path["explicit_kpoints_rel"],
        with_eigenvectors=True,
        with_group_velocities=True,
    )
    q_dict = phonon.get_qpoints_dict()
    eigenvectors = q_dict["eigenvectors"]
    group_velocities = q_dict["group_velocities"]
    if eigenvectors is None or group_velocities is None:
        raise RuntimeError("phonopy q-point run lacks eigenvectors/group velocities")

    # phonopy stores eigenvectors as (n_q, 3N, n_bands) with the column index = band and
    # the row index = 3 * atom + direction; reorder to (n_q, n_bands, N, 3) and split
    # complex values into [re, im] so the JSON matches band.yaml
    primitive = phonon.primitive
    n_atoms = len(primitive)
    n_q = len(labels)
    if eigenvectors.shape != (n_q, 3 * n_atoms, 3 * n_atoms):
        raise ValueError(
            f"eigenvector shape {eigenvectors.shape} does not match {n_q} q-points and "
            f"{n_atoms} primitive atoms"
        )
    by_band = eigenvectors.transpose(0, 2, 1).reshape(n_q, 3 * n_atoms, n_atoms, 3)
    eigenvector_array = np.stack((by_band.real, by_band.imag), axis=-1)
    norm_sq = np.sum(eigenvector_array**2, axis=(2, 3, 4))
    if not np.allclose(norm_sq, 1, atol=1e-8):
        raise RuntimeError(
            f"eigenvectors are not normalized: |e|^2 in "
            f"[{norm_sq.min():.6g}, {norm_sq.max():.6g}]"
        )

    phonon.run_mesh(list(mesh))
    phonon.run_total_dos()
    dos = phonon.get_total_dos_dict()
    phonon.run_thermal_properties(
        t_min=THERMAL_T_MIN, t_max=THERMAL_T_MAX, t_step=THERMAL_T_STEP
    )
    thermal = phonon.get_thermal_properties_dict()

    return {
        "schema_version": HARMONIC_SCHEMA_VERSION,
        "frequency_unit": FREQUENCY_UNIT,
        "primitive": {
            "lattice": np.asarray(primitive.cell),
            "symbols": list(primitive.symbols),
            "masses": np.asarray(primitive.masses),
            "frac_coords": np.asarray(primitive.scaled_positions),
        },
        "band_path": {
            # inclusive q-point index ranges with seekpath labels (GAMMA, X, SIGMA_0...)
            "segments": [
                {
                    "start_index": start,
                    "end_index": end - 1,
                    "start_label": labels[start],
                    "end_label": labels[end - 1],
                }
                for start, end in path["explicit_segments"]
            ],
            "q_points": np.asarray(path["explicit_kpoints_rel"]),
            # 1/A with 2 pi (seekpath convention; phonopy's band.yaml omits the 2 pi)
            "distances": np.asarray(path["explicit_kpoints_linearcoord"]),
            "frequencies": q_dict["frequencies"],
            "eigenvectors": np.round(eigenvector_array, EIGENVECTOR_DECIMALS),
            # THz * Angstrom
            "group_velocities": group_velocities,
        },
        "dos": {
            "mesh": [int(count) for count in mesh],
            "frequencies": dos["frequency_points"],
            # states per THz per primitive cell (phonopy tetrahedron method)
            "densities": dos["total_dos"],
        },
        "thermal_properties": {
            "temperatures": thermal["temperatures"],  # K
            "free_energy": thermal["free_energy"],  # kJ/mol
            "entropy": thermal["entropy"],  # J/K/mol
            "heat_capacity": thermal["heat_capacity"],  # J/K/mol
        },
    }


def harmonic_data_from_phono3py(ph3: Phono3py) -> dict[str, Any]:
    """Harmonic phonon data for a phono3py object with FC2 on the benchmark mesh."""
    if ph3.mesh_numbers is None:
        raise ValueError("phono3py has no mesh_numbers for the DOS calculation")
    return harmonic_phonon_data(
        phonopy_from_phono3py(ph3),
        mesh=[int(count) for count in np.asarray(ph3.mesh_numbers).ravel()],
    )
