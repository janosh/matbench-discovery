"""Perturb atomic coordinates of a pymatgen structure and analyze symmetry."""

import warnings

import numpy as np
import pandas as pd
import spglib
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2022-12-02"

rng = np.random.default_rng(0)  # ensure reproducible structure perturbations


def perturb_structure(struct: Structure, gamma: float = 1.5) -> Structure:
    """Perturb the atomic coordinates of a pymatgen structure. Used for CGCNN+P
    training set augmentation.

    Not identical but very similar to the perturbation method used in
    https://nature.com/articles/s41524-022-00891-8#Fig5.

    Args:
        struct (Structure): pymatgen structure to be perturbed
        gamma (float, optional): Weibull distribution parameter. Defaults to 1.5.

    Returns:
        Structure: Perturbed structure
    """
    perturbed = struct.copy()
    for site in perturbed:
        magnitude = rng.weibull(gamma)
        vec = rng.normal(3)  # TODO maybe make func recursive to deal with 0-vector
        vec /= np.linalg.norm(vec)  # unit vector
        site.coords += vec * magnitude
        site.to_unit_cell(in_place=True)

    return perturbed


def analyze_symmetry(structures: list[Structure]) -> pd.DataFrame:
    """Analyze symmetry of a list of structures using spglib.

    Args:
        structures (list[Structure]): List of pymatgen Structure objects

    Returns:
        pd.DataFrame: DataFrame containing symmetry information for each structure
    """
    sym_key_map = {
        "number": Key.spg_num,
        "hall_number": Key.hall_num,
        "international": MbdKey.international_spg_name,
        "hall": Key.hall_symbol,
        "choice": Key.choice_symbol,
        "pointgroup": Key.point_group,
        "wyckoffs": Key.wyckoff_symbols,
    }

    results = []
    for struct in structures:
        cell = (struct.lattice.matrix, struct.frac_coords, struct.atomic_numbers)
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=spglib.spglib.SpglibError)
            sym_data: spglib.SpglibDataset = spglib.get_symmetry_dataset(cell)

        sym_info = {
            new_key: getattr(sym_data, old_key)
            for old_key, new_key in sym_key_map.items()
        }
        sym_info[Key.n_rot_ops] = len(sym_data.rotations)
        sym_info[Key.n_trans_ops] = len(sym_data.translations)
        sym_info[Key.n_sym_ops] = sym_info[Key.n_rot_ops] + sym_info[Key.n_trans_ops]

        results.append(sym_info)

    df_sym = pd.DataFrame(results)
    df_sym.index.name = Key.mat_id
    return df_sym


def compare_symmetry(
    df_ml: pd.DataFrame,
    df_dft: pd.DataFrame,
    ml_structs: dict[str, Structure],
    dft_structs: dict[str, Structure],
) -> pd.DataFrame:
    """Compare symmetry information between ML model and DFT reference data.

    Args:
        df_ml (pd.DataFrame): symmetry information for ML model
        df_dft (pd.DataFrame): symmetry information for DFT reference
        ml_structs (dict[str, Structure]): Map material IDs to ML structures
        dft_structs (dict[str, Structure]): Map material IDs to DFT structures

    Returns:
        pd.DataFrame: with added columns for symmetry differences
    """
    df_combined = df_ml.copy()

    # Calculate differences
    df_combined[MbdKey.spg_num_diff] = df_combined[Key.spg_num] - df_dft[Key.spg_num]
    df_combined[MbdKey.n_sym_ops_diff] = (
        df_combined[Key.n_sym_ops] - df_dft[Key.n_sym_ops]
    )

    structure_matcher = StructureMatcher()
    common_ids = set(ml_structs) & set(dft_structs)
    df_combined[[Key.rmsd, Key.max_pair_dist]] = [
        structure_matcher.get_rms_dist(ml_structs[mat_id], dft_structs[mat_id])
        for mat_id in common_ids
    ]

    return df_combined


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    gamma = 1.5
    samples = np.array([rng.weibull(gamma) for _ in range(10_000)])
    mean = samples.mean()

    # reproduces the dist in https://www.nature.com/articles/s41524-022-00891-8#Fig5
    ax = plt.hist(samples, bins=100)
    # add vertical line at the mean
    plt.axvline(mean, color="gray", linestyle="dashed", linewidth=1)
    # annotate the mean line
    plt.annotate(
        f"{mean=:.2f}",
        xy=(mean, 1),
        # use ax coords for y
        xycoords=("data", "axes fraction"),
        # add text offset
        xytext=(10, -20),
        textcoords="offset points",
    )
