"""This package contains phonon-related functionality."""

from typing import Final, Literal, NotRequired, TypedDict

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.phonons.calc_kappa import calc_kappa_for_structure


def read_kappa_json(path: str) -> pd.DataFrame:
    """Read a kappa prediction/reference JSON file indexed by material_id.

    Renames the ID column for submissions that use mp_id instead of material_id.
    """
    df_kappa = pd.read_json(path)
    columns = list(df_kappa)
    if "mp_id" in columns:
        df_kappa = df_kappa.rename(columns={"mp_id": Key.mat_id})
    elif Key.mat_id not in columns:
        raise ValueError(
            f"read_kappa_json requires an 'mp_id' or {Key.mat_id!r} column, "
            f"got {columns=}"
        )
    return df_kappa.set_index(Key.mat_id)


class KappaCalcParams(TypedDict):
    """Parameters for thermal conductivity calculation across all models.

    This TypedDict provides type safety for parameters passed to kappa calculation
    functions like calc_kappa_for_structure.

    Required parameters (all models):
        displacement_distance: Displacement distance for phono3py finite differences (Å)
        temperatures: Temperatures in Kelvin for conductivity calculation
        ase_optimizer: ASE optimizer name (e.g., 'FIRE', 'BFGS', 'LBFGS')
        max_steps: Maximum relaxation steps
        force_max: Maximum force convergence criterion (eV/Å)
        symprec: Symmetry precision for spglib
        enforce_relax_symm: Whether to enforce symmetry during relaxation
        save_forces: Whether to save force sets to disk
        out_dir: Output directory for results

    Optional parameters (model-specific):
        ase_filter: Cell filter for relaxation ('frechet' or 'exp') - NequIP, Allegro
        checkpoint: Model checkpoint path or URL - MACE
        conductivity_broken_symm: Calculate kappa if symmetry breaks - MACE
    """

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
    # Return True if all frequencies are NaN, indicating invalid or missing data
    if np.all(pd.isna(frequencies)):
        return True

    # Check for imaginary frequencies in non-acoustic modes at gamma point (q=0)
    # Indices 3+ correspond to optical modes which should never be negative
    if np.any(frequencies[0, 3:] < 0):
        return True

    # Check acoustic modes at gamma point against threshold. First 3 modes at q=0
    # are acoustic and may be slightly negative due to numerical noise
    if np.any(frequencies[0, :3] < threshold):
        return True

    # Check for imaginary frequencies at any q-point except gamma
    # All frequencies should be positive away from gamma point
    return bool(np.any(frequencies[1:] < 0))
