"""Shared helpers for repository script and pipeline tests."""

import importlib.util
import sys
from types import ModuleType
from typing import Any

import numpy as np

from matbench_discovery import ROOT


def import_repo_script(module_name: str, rel_path: str) -> ModuleType:
    """Import a repository-local script without package-name collisions."""
    spec = importlib.util.spec_from_file_location(module_name, f"{ROOT}/{rel_path}")
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {rel_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    return module


def make_harmonic_record(
    n_atoms: int = 2, n_q: int = 4, material_id: str = "material-1"
) -> dict[str, Any]:
    """Build a minimal harmonic phonon sidecar record with a two-segment band path.

    Mirrors ``matbench_discovery.phonons.harmonic.harmonic_phonon_data`` output:
    normalized eigenvectors at every q-point, segments GAMMA-X and X-L splitting the
    path in half, thermal properties at 0 and 300 K.
    """
    n_bands = 3 * n_atoms
    eigenvectors = np.zeros((n_q, n_bands, n_atoms, 3, 2))
    eigenvectors[..., 0] = 1 / np.sqrt(3 * n_atoms)
    half = n_q // 2
    return {
        "material_id": material_id,
        "schema_version": 1,
        "frequency_unit": "THz",
        "errors": {},
        "primitive": {
            "lattice": np.eye(3) * 3.123456789,
            "symbols": ["Na", "Cl", "K", "Br"][:n_atoms],
            "masses": [22.99, 35.45, 39.1, 79.9][:n_atoms],
            "frac_coords": np.linspace(0, 0.75, n_atoms)[:, None] * np.ones(3),
        },
        "band_path": {
            "segments": [
                {
                    "start_index": 0,
                    "end_index": half - 1,
                    "start_label": "GAMMA",
                    "end_label": "X",
                },
                {
                    "start_index": half,
                    "end_index": n_q - 1,
                    "start_label": "X",
                    "end_label": "L",
                },
            ],
            "q_points": np.linspace(0, 0.5, n_q)[:, None] * np.ones(3),
            "distances": np.linspace(0, 1.23456789, n_q),
            "frequencies": np.arange(n_q * n_bands).reshape(n_q, n_bands) / 7,
            "eigenvectors": eigenvectors,
            "group_velocities": np.zeros((n_q, n_bands, 3)),
        },
        "dos": {"mesh": [2, 2, 2], "frequencies": [0.0, 1.0], "densities": [0, 1]},
        "thermal_properties": {
            "temperatures": [0.0, 300.0],
            "free_energy": [1.0, 0.0],
            "entropy": [0.0, 1.0],
            "heat_capacity": [0.0, 1.0],
        },
    }
