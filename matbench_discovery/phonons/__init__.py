"""This package contains phonon-related functionality."""

from typing import TYPE_CHECKING, Literal, NotRequired, TypedDict

import numpy as np
import pandas as pd
from pymatviz.enums import Key

if TYPE_CHECKING:
    from matbench_discovery.phonons.thermal_conductivity import (
        calc_kappa_for_structure as calc_kappa_for_structure,
    )


def __getattr__(name: str) -> object:
    """Lazily expose phono3py-backed helpers without requiring the `phonons` extra."""
    if name == "calc_kappa_for_structure":
        from matbench_discovery.phonons.thermal_conductivity import (
            calc_kappa_for_structure,
        )

        return calc_kappa_for_structure
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def read_kappa_json(path: str) -> pd.DataFrame:
    """Read a kappa prediction/reference JSON file indexed by material_id.

    Renames the ID column for submissions that use mp_id instead of material_id.
    """
    df_kappa = pd.read_json(path)
    if "mp_id" in df_kappa:
        df_kappa = df_kappa.rename(columns={"mp_id": Key.mat_id})
    if Key.mat_id in df_kappa:
        return df_kappa.set_index(Key.mat_id)
    raise ValueError(
        f"read_kappa_json requires an 'mp_id' or {Key.mat_id!r} column, "
        f"got columns={list(df_kappa)!r}"
    )


class KappaCalcParams(TypedDict):
    """Shared keyword arguments for phonon thermal conductivity calculations."""

    # Required for all models
    displacement_distance: float
    temperatures: list[float]
    ase_optimizer: str
    max_steps: int
    force_max: float
    symprec: float
    enforce_relax_symm: bool
    save_forces: bool
    out_dir: str
    # Optional model-specific parameters
    ase_filter: NotRequired[Literal["frechet", "exp"]]
    checkpoint: NotRequired[str]
    conductivity_broken_symm: NotRequired[bool]


def check_imaginary_freqs(frequencies: np.ndarray, threshold: float = -0.01) -> bool:
    """Check if frequencies are imaginary.

    Args:
        frequencies (np.ndarray): Frequencies to check
        threshold (float): Threshold for imaginary frequencies. Defaults to -0.01.

    Returns:
        bool: True if imaginary frequencies are found.
    """
    return bool(
        # Return True if all frequencies are NaN, indicating invalid or missing data
        np.all(pd.isna(frequencies))
        # Check for imaginary frequencies in non-acoustic modes at gamma point (q=0)
        # Indices 3+ correspond to optical modes which should never be negative
        or np.any(frequencies[0, 3:] < 0)
        # Check acoustic modes at gamma point against threshold. First 3 modes at q=0
        # are acoustic and may be slightly negative due to numerical noise
        or np.any(frequencies[0, :3] < threshold)
        # Check for imaginary frequencies at any q-point except gamma
        # All frequencies should be positive away from gamma point
        or np.any(frequencies[1:] < 0)
    )
