"""This package contains phonon-related functionality."""

import numpy as np
import pandas as pd
from pymatviz.enums import Key


def read_kappa_json(path: str) -> pd.DataFrame:
    """Read a kappa prediction/reference JSON file indexed by material_id.

    Requires canonical field names and unique, non-empty material IDs.
    """
    from matbench_discovery.phonons.schema import normalize_kappa_dataframe

    df_kappa = pd.read_json(path)
    if str(Key.mat_id) not in df_kappa:
        raise ValueError(
            "read_kappa_json requires a canonical 'material_id' column, "
            f"got columns={list(df_kappa)!r}"
        )
    df_kappa = normalize_kappa_dataframe(df_kappa)
    return df_kappa.set_index(str(Key.mat_id), verify_integrity=True)


# acoustic modes at gamma come out slightly negative from numerical noise (measured
# ~-5e-8 THz on stable fcc Cu), so only frequencies below this count as imaginary
IMAGINARY_FREQ_THRESHOLD = -0.01


def check_imaginary_freqs(
    frequencies: np.ndarray, threshold: float = IMAGINARY_FREQ_THRESHOLD
) -> bool:
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
