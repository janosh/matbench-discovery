"""Harmonic phonon data derived from second-order force constants.

The kappa pipeline only needs mesh frequencies, but the same FC2 also gives the full
harmonic picture: high-symmetry band path with eigenvectors and group velocities,
tetrahedron DOS on the benchmark mesh, and thermal properties. The sidecar retains
unrounded FC2 and its defining cell data so new metrics can be evaluated offline.

Eigenvector layout follows phonopy's ``band.yaml`` convention: dynamical-matrix
eigenvectors normalized to ``sum |e|^2 = 1``, as ``[re, im]`` pairs per atom and
Cartesian direction (the site animates them with matterviz ``PhononModeExplorer``).
They are kept at every path q-point, not just the labeled ones the site animates
(tens of MB gzipped per model, on the order of the prediction file), so bands can
later be tracked through crossings by eigenvector overlap.

For a new metric, load a record from ``*-phonons.json.gz`` and reconstruct it with
``phonons_from_harmonic_data(record)`` (a ``matscikit.Phonons``). Evaluate the required
q-points or mesh on that object, then keep scoring functions in
``metrics/phonons.py``. Existing mesh spectra already carry q-point multiplicities;
band samples are not a Brillouin-zone measure. Use the FC2 input for
precision-sensitive work, since displayed eigenvectors are rounded. New derived
properties belong in separate analysis functions, with their own entry in
``errors``; they must never enter the conductivity ``errors`` list.
The run-info sidecar supplies the source, checkpoint, and package provenance. Old
v1 artifacts lack FC2, v2 artifacts store FC2 in phonopy's supercell atom order and
v3 full FC2; none can be replayed, so regenerate them explicitly. v4 stores compact
FC2 (primitive-atom rows) with the supercell fractional coordinates, making the FC2
row order self-describing.
"""

from __future__ import annotations

import traceback
from collections.abc import Callable, Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
from matscikit import Phonons, Structure, band_path

from matbench_discovery.phonons import IMAGINARY_FREQ_THRESHOLD

if TYPE_CHECKING:
    from collections.abc import Sequence

    from matscikit import MeshResult

# v4: compact FC2 (primitive-atom rows) in matscikit's recorded supercell atom order,
# v3: full FC2 in that order, v2: full FC2 in phonopy's order
HARMONIC_SCHEMA_VERSION = 4
# HPKOT q-point spacing along the path in 1/A (2 pi convention): the FCC Cu
# GAMMA-X segment gets 69 points, in the range of phonopy's 51-point band.conf default
BAND_PATH_REFERENCE_DISTANCE = 0.025
FREQUENCY_UNIT = "THz"
THERMAL_TEMPERATURES = np.arange(0.0, 1001.0, 10.0)  # K
# eigenvector components are rounded to this many decimals: matterviz re-normalizes
# within 1e-3 of |e|^2 = 1 and the rounding keeps the sidecar a few MB per material
EIGENVECTOR_DECIMALS = 5


def harmonic_band_data(phonons: Phonons) -> dict[str, Any]:
    """HPKOT band path with frequencies, eigenvectors and group velocities.

    The path is computed on the primitive cell in its own basis (q-points fractional
    in the primitive reciprocal basis); labels follow HPKOT (GAMMA, X, SIGMA_0, ...).
    """
    primitive = phonons.primitive
    path = band_path(
        primitive,
        symprec=phonons.symprec,
        reference_distance=BAND_PATH_REFERENCE_DISTANCE,
    )
    labels = path.labels
    modes = phonons.modes(path.q_points, with_group_velocities=True)
    eigenvectors = modes["eigenvectors"]
    n_atoms = len(primitive)
    n_q = len(labels)
    if eigenvectors.shape != (n_q, 3 * n_atoms, 3 * n_atoms):
        raise ValueError(
            f"eigenvector shape {eigenvectors.shape} does not match {n_q} q-points and "
            f"{n_atoms} primitive atoms"
        )
    # eigenvectors are (n_q, 3N, n_bands) with column = band and row = 3 * atom +
    # direction; reorder to (n_q, n_bands, N, 3) and split complex values into
    # [re, im] so the JSON matches band.yaml
    by_band = eigenvectors.transpose(0, 2, 1).reshape(n_q, 3 * n_atoms, n_atoms, 3)
    eigenvector_array = np.stack((by_band.real, by_band.imag), axis=-1)
    norm_sq = np.sum(eigenvector_array**2, axis=(2, 3, 4))
    # rtol=0: allclose's default rtol=1e-5 against a desired value of 1 would loosen
    # this to 1e-5; real eigenvectors are normalized to ~1e-15
    if not np.allclose(norm_sq, 1, atol=1e-8, rtol=0):
        raise RuntimeError(
            f"eigenvectors are not normalized: |e|^2 in "
            f"[{norm_sq.min():.6g}, {norm_sq.max():.6g}]"
        )

    return {
        # inclusive q-point index ranges with HPKOT labels (GAMMA, X, SIGMA_0...)
        "segments": [
            {
                "start_index": start,
                "end_index": end - 1,
                "start_label": labels[start],
                "end_label": labels[end - 1],
            }
            for start, end in path.segments
        ],
        "q_points": path.q_points,
        # 1/A with 2 pi (HPKOT convention; phonopy's band.yaml omits the 2 pi)
        "distances": path.distances,
        "frequencies": modes["frequencies"],
        "eigenvectors": np.round(eigenvector_array, EIGENVECTOR_DECIMALS),
        # THz * Angstrom
        "group_velocities": modes["group_velocities"],
    }


def harmonic_thermal_data(mesh_result: MeshResult) -> dict[str, Any]:
    """Thermal curves of a ``Phonons.mesh(..., temperatures_k=...)`` result, with the
    count of imaginary modes excluded.

    Modes below 1e-3 THz are excluded (imaginary modes and the acoustic modes at
    Gamma), so for a dynamically unstable structure every curve below is computed
    from a subset of the modes. Record the count so consumers can say so rather than
    presenting the truncated values as the material's thermodynamics. Counted
    against the same threshold as check_imaginary_freqs: modes from the -0.01..0 THz
    noise band are physically zero and would flag stable materials.
    """
    if mesh_result.temperatures_k is None:
        raise ValueError(
            "harmonic_thermal_data needs a mesh result computed with temperatures_k"
        )
    mesh_frequencies = mesh_result.frequencies
    return {
        # modes excluded from every series below (0 = all modes contribute)
        "n_imaginary_modes_excluded": int(
            (mesh_frequencies < IMAGINARY_FREQ_THRESHOLD).sum()
        ),
        "n_mesh_modes": int(mesh_frequencies.size),
        "temperatures": mesh_result.temperatures_k,
        "free_energy": mesh_result.free_energy_kj_per_mol,
        "entropy": mesh_result.entropy_j_per_k_mol,
        "heat_capacity": mesh_result.heat_capacity_j_per_k_mol,
    }


def phonons_from_harmonic_data(data: Mapping[str, Any]) -> Phonons:
    """Restore benchmark FC2 for new metrics without a model or force calls.

    The artifact stores unrounded compact force constants in eV/Å² (primitive-atom
    rows, expanded to the full supercell by its lattice translations) with the cell
    data that defines them; the recorded supercell coordinates must match the rebuilt
    supercell row for row. No force-constant fitting or symmetrization is performed.
    """
    if data["schema_version"] != HARMONIC_SCHEMA_VERSION:
        raise ValueError(f"Unsupported harmonic schema {data['schema_version']!r}")
    state = data["fc2"]
    unit = state["unitcell"]
    unitcell = Structure(
        np.asarray(unit["cell"]),
        list(unit["symbols"]),
        np.asarray(unit["scaled_positions"]),
        fractional=True,
    )
    phonons = Phonons(
        unitcell,
        np.asarray(state["supercell_matrix"]),
        primitive_matrix=np.asarray(state["primitive_matrix"]),
        symprec=float(state["symprec"]),
    )
    stored = np.asarray(state["supercell_frac_coords"])
    diff = phonons.supercell.frac_coords - stored
    if diff.shape != stored.shape or np.abs(diff - np.round(diff)).max() > 1e-6:
        raise ValueError(
            "stored FC2 supercell does not match the rebuilt supercell order"
        )
    phonons.set_force_constants(np.asarray(state["force_constants"], dtype=float))
    return phonons


def harmonic_phonon_data(phonons: Phonons, *, mesh: Sequence[int]) -> dict[str, Any]:
    """Save reusable benchmark FC2 and independently evaluate harmonic properties.

    Add new properties as ordinary analysis functions using the restored Phonons
    object. Each analysis owns its errors; failed band paths cannot discard FC2,
    mesh spectra, DOS, or thermodynamics. DOS and thermodynamics share the mesh's
    diagonalization, so they fail with it. The pipeline records these failures
    separately from kappa.
    """
    # full FC2 repeats the primitive rows under lattice translations (~1.4 GB/model)
    force_constants = phonons.compact_force_constants
    if not np.isfinite(force_constants).all():
        raise ValueError("Harmonic data requires finite force constants")
    primitive, unitcell = phonons.primitive, phonons.unitcell
    data: dict[str, Any] = {
        "schema_version": HARMONIC_SCHEMA_VERSION,
        "frequency_unit": FREQUENCY_UNIT,
        "fc2": {
            "unitcell": {
                "cell": np.asarray(unitcell.lattice),
                "symbols": list(unitcell.symbols),
                "masses": np.asarray(unitcell.masses),
                "scaled_positions": np.asarray(unitcell.frac_coords),
            },
            "supercell_matrix": np.asarray(phonons.supercell_matrix),
            "primitive_matrix": np.asarray(phonons.primitive_matrix),
            "symprec": phonons.symprec,
            "supercell_frac_coords": np.asarray(phonons.supercell.frac_coords),
            "force_constants": force_constants,
        },
        "primitive": {
            "lattice": np.asarray(primitive.lattice),
            "symbols": list(primitive.symbols),
            "masses": np.asarray(primitive.masses),
            "frac_coords": np.asarray(primitive.frac_coords),
        },
        "errors": {},
    }

    def collect(name: str, calculate: Callable[[], dict[str, Any]]) -> None:
        """Store the sections ``calculate`` returns, or its failure with traceback."""
        try:
            data.update(calculate())
        except Exception as exc:  # noqa: BLE001 - preserve other analyses and FC2
            data["errors"][name] = {
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }

    collect("band_path", lambda: {"band_path": harmonic_band_data(phonons)})
    collect("mesh", lambda: harmonic_mesh_data(phonons, mesh))
    return data


def harmonic_mesh_data(phonons: Phonons, mesh: Sequence[int]) -> dict[str, Any]:
    """Mesh frequencies with multiplicities, tetrahedron DOS (states per THz per
    primitive cell) and thermal properties from one diagonalization.
    """
    result = phonons.mesh(
        list(mesh), total_dos=True, temperatures_k=THERMAL_TEMPERATURES.tolist()
    )
    return {
        "mesh": {
            "numbers": list(mesh),
            "q_points": result.q_points,
            "weights": result.weights,
            "frequencies": result.frequencies,
        },
        "dos": {
            "mesh": list(mesh),
            "frequencies": result.frequency_points,
            "densities": result.total_dos,
        },
        "thermal_properties": harmonic_thermal_data(result),
    }
