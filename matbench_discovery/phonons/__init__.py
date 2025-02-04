"""This package contains phonon-related functionality."""

from typing import Final

import numpy as np
import pandas as pd

# q-point mesh (which phonon modes to sample) based on international space group number
spg_num_q_mesh_map: Final[dict[int, tuple[int, int, int]]] = {
    225: (19, 19, 19),  # rocksalt
    186: (19, 19, 15),  # wurtzite
    216: (19, 19, 19),  # zincblende
}


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
