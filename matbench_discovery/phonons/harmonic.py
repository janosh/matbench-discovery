"""Harmonic phonon data derived from second-order force constants.

The kappa pipeline only needs mesh frequencies, but the same FC2 also gives the full
harmonic picture: high-symmetry band path with eigenvectors and group velocities,
tetrahedron DOS on the benchmark mesh, and thermal properties. The sidecar retains
unrounded FC2 and its defining cell data so new metrics can be evaluated offline.

Eigenvector layout follows phonopy's ``band.yaml``: mass-weighted dynamical-matrix
eigenvectors normalized to ``sum |e|^2 = 1``, as ``[re, im]`` pairs per atom and
Cartesian direction (the site animates them with matterviz ``PhononModeExplorer``).
They are kept at every path q-point, not just the labeled ones the site animates
(tens of MB gzipped per model, on the order of the prediction file), so bands can
later be tracked through crossings by eigenvector overlap.

For a new metric, load a record from ``*-phonons.json.gz`` and reconstruct it with
``phonopy_from_harmonic_data(record)``. Evaluate the required q-points or mesh on that
object, then keep scoring functions in ``metrics/phonons.py``. Existing mesh spectra
already carry q-point multiplicities; band samples are not a Brillouin-zone measure.
Use the FC2 input for precision-sensitive work, since displayed eigenvectors are
rounded. New derived properties belong in separate analysis functions, with their
own entry in ``errors``; they must never enter the conductivity ``errors`` list.
The run-info sidecar supplies the source, checkpoint, and package provenance. Old
v1 harmonic artifacts lack FC2 and cannot be replayed; regenerate them explicitly.
"""

from __future__ import annotations

import traceback
from collections.abc import Callable, Mapping
from typing import TYPE_CHECKING, Any

import numpy as np

from matbench_discovery.phonons import IMAGINARY_FREQ_THRESHOLD

if TYPE_CHECKING:
    from collections.abc import Sequence

    from phono3py.api_phono3py import Phono3py
    from phonopy import Phonopy

HARMONIC_SCHEMA_VERSION = 2
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

    # reuse phonopy's own tolerance: at seekpath's stricter 1e-5 default a run
    # configured with a looser symprec can standardize to a different primitive cell,
    # which trips the unimodular basis-change guard below
    path = seekpath.get_explicit_k_path(
        phonon.primitive.totuple(),
        reference_distance=BAND_PATH_REFERENCE_DISTANCE,
        symprec=phonon.primitive_symmetry.tolerance,
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


def harmonic_band_data(phonon: Phonopy) -> dict[str, Any]:
    """Compute seekpath bands, mass-weighted eigenvectors, and group velocities."""
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
    # rtol=0: allclose's default rtol=1e-5 against a desired value of 1 would loosen
    # this to 1e-5. Measured |norm^2 - 1| is 1.8e-15 (8x f64 eps) on fcc Cu/Al/Au, so
    # the 1e-8 absolute bound still clears real eigenvectors by ~6 orders of magnitude
    if not np.allclose(norm_sq, 1, atol=1e-8, rtol=0):
        raise RuntimeError(
            f"eigenvectors are not normalized: |e|^2 in "
            f"[{norm_sq.min():.6g}, {norm_sq.max():.6g}]"
        )

    return {
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
    }


def harmonic_thermal_data(phonon: Phonopy) -> dict[str, Any]:
    """Compute thermal curves on an already evaluated mesh, with instability counts."""
    phonon.run_thermal_properties(
        t_min=THERMAL_T_MIN, t_max=THERMAL_T_MAX, t_step=THERMAL_T_STEP
    )
    thermal = phonon.get_thermal_properties_dict()
    # phonopy's ThermalProperties silently drops modes below cutoff_frequency=0, so for
    # a dynamically unstable structure every curve below is computed from a subset of
    # the modes (measured: Cv(1000 K) = 17.7 vs the 49.9 J/K/mol classical limit on a
    # 2-atom cell with imaginary branches). Record the count so consumers can say so
    # rather than presenting the truncated values as the material's thermodynamics.
    # Counted against the same threshold as check_imaginary_freqs: modes dropped from
    # the -0.01..0 THz noise band are physically zero and would flag stable materials
    mesh_frequencies = np.asarray(phonon.get_mesh_dict()["frequencies"])
    n_imaginary_modes = int((mesh_frequencies < IMAGINARY_FREQ_THRESHOLD).sum())

    return {
        # modes excluded from every series below (0 = all modes contribute)
        "n_imaginary_modes_excluded": n_imaginary_modes,
        "n_mesh_modes": int(mesh_frequencies.size),
        "temperatures": thermal["temperatures"],  # K
        "free_energy": thermal["free_energy"],  # kJ/mol
        "entropy": thermal["entropy"],  # J/K/mol
        "heat_capacity": thermal["heat_capacity"],  # J/K/mol
    }


def phonopy_from_harmonic_data(data: Mapping[str, Any]) -> Phonopy:
    """Restore benchmark FC2 for new metrics without a calculator or force calls.

    The artifact stores unrounded force constants in eV/Å² and explicit cell bases,
    masses, symmetry tolerance, frequency conversion, and optional NAC parameters.
    No force-constant fitting or symmetrization is performed during reconstruction.
    """
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    if data["schema_version"] != HARMONIC_SCHEMA_VERSION:
        raise ValueError(f"Unsupported harmonic schema {data['schema_version']!r}")
    state = data["fc2"]
    phonon = Phonopy(
        PhonopyAtoms(**state["unitcell"]),
        supercell_matrix=state["supercell_matrix"],
        primitive_matrix=state["primitive_matrix"],
        symprec=state["symprec"],
    )
    phonon.unit_conversion_factor = state["frequency_factor"]
    phonon.force_constants = np.asarray(state["force_constants"], dtype=float)
    phonon.nac_params = state["nac_params"]
    return phonon


def harmonic_phonon_data(phonon: Phonopy, *, mesh: Sequence[int]) -> dict[str, Any]:
    """Save reusable benchmark FC2 and independently evaluate harmonic properties.

    Add new properties as ordinary analysis functions using the restored Phonopy
    object. Each analysis owns its errors; failed band paths cannot discard FC2,
    mesh spectra, DOS, or thermodynamics. Only mesh-dependent analyses stop when
    the mesh fails. The pipeline records these failures separately from kappa.
    """
    if phonon.force_constants is None or not np.isfinite(phonon.force_constants).all():
        raise ValueError("Harmonic data requires finite force constants")
    primitive, unitcell = phonon.primitive, phonon.unitcell
    data: dict[str, Any] = {
        "schema_version": HARMONIC_SCHEMA_VERSION,
        "frequency_unit": FREQUENCY_UNIT,
        "fc2": {
            "unitcell": {
                "cell": np.asarray(unitcell.cell),
                "symbols": list(unitcell.symbols),
                "masses": np.asarray(unitcell.masses),
                "scaled_positions": np.asarray(unitcell.scaled_positions),
                "magnetic_moments": unitcell.magnetic_moments,
            },
            "supercell_matrix": phonon.supercell_matrix,
            "primitive_matrix": phonon.primitive_matrix,
            "symprec": phonon.symmetry.tolerance,
            "frequency_factor": phonon.unit_conversion_factor,
            "force_constants": phonon.force_constants,
            "nac_params": phonon.nac_params,
        },
        "primitive": {
            "lattice": np.asarray(primitive.cell),
            "symbols": list(primitive.symbols),
            "masses": np.asarray(primitive.masses),
            "frac_coords": np.asarray(primitive.scaled_positions),
        },
        "errors": {},
    }

    def collect(name: str, calculate: Callable[[], dict[str, Any]]) -> None:
        """Retain independent analysis failures with their traceback for diagnosis."""
        try:
            data[name] = calculate()
        except Exception as exc:  # noqa: BLE001 - preserve other analyses and FC2
            data["errors"][name] = {
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }

    collect("band_path", lambda: harmonic_band_data(phonon))

    def mesh_data() -> dict[str, Any]:
        """Store q-points and their multiplicities together with mesh frequencies."""
        phonon.run_mesh(list(mesh))
        mesh_result = phonon.get_mesh_dict()
        return {
            "numbers": list(mesh),
            "q_points": mesh_result["qpoints"],
            "weights": mesh_result["weights"],
            "frequencies": mesh_result["frequencies"],
        }

    collect("mesh", mesh_data)
    if "mesh" not in data:
        return data

    def dos_data() -> dict[str, Any]:
        """Compute tetrahedron DOS in states per THz per primitive cell."""
        phonon.run_total_dos()
        dos = phonon.get_total_dos_dict()
        return {
            "mesh": list(mesh),
            "frequencies": dos["frequency_points"],
            "densities": dos["total_dos"],
        }

    collect("dos", dos_data)
    collect("thermal_properties", lambda: harmonic_thermal_data(phonon))
    return data


def harmonic_data_from_phono3py(ph3: Phono3py) -> dict[str, Any]:
    """Harmonic phonon data for a phono3py object with FC2 on the benchmark mesh."""
    if ph3.mesh_numbers is None:
        raise ValueError("phono3py has no mesh_numbers for the DOS calculation")
    return harmonic_phonon_data(
        phonopy_from_phono3py(ph3),
        mesh=[int(count) for count in np.asarray(ph3.mesh_numbers).ravel()],
    )
