"""Calculate MLFF pair-repulsion curves for diatomic molecules.

Thanks to Tamas Stenczel who first did this type of PES smoothness and physicality
analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP paper
https://arxiv.org/abs/2401.00096 (see fig. 56).
"""

# %%
from __future__ import annotations

from typing import TYPE_CHECKING, Literal, get_args

import numpy as np
from ase import Atoms
from ase.data import chemical_symbols
from tqdm import tqdm

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase.calculators.calculator import Calculator

__date__ = "2024-03-31"
DiatomicsType = Literal["homo-nuclear", "hetero-nuclear"]
homo_nuc = get_args(DiatomicsType)[0]
atom_num_symbol_map = dict(enumerate(chemical_symbols, start=0))
# one curve's distances/energies/forces lists, keyed by formula in DiatomicResults
CurveDict = dict[str, list[float | list[list[float]]]]
DiatomicResults = dict[str, CurveDict]


def generate_diatomics(
    elem1: str,
    elem2: str,
    distances: Sequence[float | int] | np.ndarray,
    box_size: float = 50.0,
) -> list[Atoms]:
    """Build diatomic molecules centered in a large cubic cell for given distances.

    A large periodic box (rather than ``pbc=False`` with no cell) is used so calculators
    that require an invertible cell matrix (graph MLIPs like CHGNet, eqnorm, TACE, which
    otherwise raise "box matrix is not invertible" / "Singular matrix") work. The box is
    large enough that periodic images don't interact within any short-range MLIP cutoff,
    so energies/forces match the isolated-molecule values.

    Args:
        elem1 (str): Chemical symbol of the first element.
        elem2 (str): Chemical symbol of the second element.
        distances (Sequence[float]): Distances to sample at.
        box_size (float): Cubic cell edge length in Å. Defaults to 50.

    Returns:
        list[Atoms]: Diatomic molecules with elements at given distances.
    """
    center = box_size / 2
    if len(distances) and (max_dist := max(distances)) >= center:
        raise ValueError(
            f"{max_dist=:.3f} Å must be < box_size/2 = {center=:.3f} Å, else "
            "periodic images give a wrong min-image distance; raise box_size"
        )
    return [
        Atoms(
            f"{elem1}{elem2}",
            positions=[[center, center, center], [center + dist, center, center]],
            cell=[box_size] * 3,
            pbc=True,
        )
        for dist in distances
    ]


def calc_diatomic_curve(
    pairs: Sequence[tuple[str | int, str | int]],
    calculator: Calculator,
    model_name: str,
    distances: Sequence[float | int] | np.ndarray,
    results: DiatomicResults,
) -> DiatomicResults:
    """Generate potential energy and forces data for diatomic molecules.

    Args:
        pairs (list[tuple[str | int, str | int]]): List of element pairs to calculate.
            Each pair can be specified as element symbols or atomic numbers.
        calculator (Calculator): ASE calculator instance.
        model_name (str): Name of the model for the output file.
        distances (list[float]): Distances to calculate potential energy at.
        results (DiatomicResults): Results dict to collect energies and forces at given
            distances for all diatomic curves. Will be updated in-place.

    Returns:
        DiatomicResults: Potential energy and forces data for diatomic molecules.
    """
    # saving results in dict: {"symbol-symbol": {"energies": [...], "forces": [...]}}
    expected_distances = np.asarray(distances)
    for idx, (z1, z2) in (pbar := tqdm(enumerate(pairs, start=1))):
        # Convert atomic numbers to symbols if needed
        elem1 = atom_num_symbol_map.get(z1, z1)
        elem2 = atom_num_symbol_map.get(z2, z2)
        formula = f"{elem1}-{elem2}"
        curve_data = results.get(formula, {})
        len_energies = len(curve_data.get("energies", []))
        len_forces = len(curve_data.get("forces", []))
        cached_distances = np.asarray(curve_data.get("distances", []))

        # skip only when cached samples match the current distance grid exactly
        if len_energies == len_forces == len(expected_distances) and np.array_equal(
            cached_distances, expected_distances
        ):
            continue

        pbar.set_description(
            f"{idx}/{len(pairs)} {formula} diatomic curve with {model_name}"
        )

        # reset curve_data in case we had prior results
        curve_data = results[formula] = {
            "distances": expected_distances.tolist(),
            "energies": [],
            "forces": [],
        }
        try:
            for atoms in generate_diatomics(elem1, elem2, distances):
                curve_data["energies"].append(calculator.get_potential_energy(atoms))
                curve_data["forces"].append(calculator.get_forces(atoms).tolist())
        except (
            ValueError,
            RuntimeError,
            KeyError,
            NotImplementedError,
            AssertionError,
        ) as exc:
            # elements outside a model's training set (EMT's supported metals, GRACE's
            # asserted element list, ...) raise here; skip that pair's curve rather than
            # aborting the whole sweep
            print(f"{idx}/{len(pairs)} {formula} failed: {exc!r}")
            results[formula] = {
                "distances": expected_distances.tolist(),
                "energies": [],
                "forces": [],
            }

    return results
