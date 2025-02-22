"""Calculate MLFF pair-repulsion curves for diatomic molecules.

Thanks to Tamas Stenczel who first did this type of PES smoothness and physicality
analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP paper
https://arxiv.org/abs/2401.00096 (see fig. 56).
"""

# %%
from __future__ import annotations

import os
from typing import TYPE_CHECKING, Literal, get_args

from ase import Atoms
from ase.data import chemical_symbols
from tqdm import tqdm

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase.calculators.calculator import Calculator

__date__ = "2024-03-31"
module_dir = os.path.dirname(__file__)
DiatomicsType = Literal["homo-nuclear", "hetero-nuclear"]
homo_nuc, hetero_nuc = get_args(DiatomicsType)
atom_num_symbol_map = dict(enumerate(chemical_symbols, start=0))


def generate_diatomics(
    elem1: str, elem2: str, distances: Sequence[float]
) -> list[Atoms]:
    """Build diatomic molecules in vacuum for given distances.

    Args:
        elem1 (str): Chemical symbol of the first element.
        elem2 (str): Chemical symbol of the second element.
        distances (Sequence[float]): Distances to sample at.

    Returns:
        list[Atoms]: Diatomic molecules with elements at given distances.
    """
    return [
        Atoms(f"{elem1}{elem2}", positions=[[0, 0, 0], [dist, 0, 0]], pbc=False)
        for dist in distances
    ]


def calc_diatomic_curve(
    pairs: Sequence[tuple[str | int, str | int]],
    calculator: Calculator,
    model_name: str,
    distances: Sequence[float],
    results: dict[str, dict[str, list[float | list[list[float]]]]],
) -> dict[str, dict[str, list[float | list[list[float]]]]]:
    """Generate potential energy and forces data for diatomic molecules.

    Args:
        pairs (list[tuple[str | int, str | int]]): List of element pairs to calculate.
            Each pair can be specified as element symbols or atomic numbers.
        calculator (Calculator): ASE calculator instance.
        model_name (str): Name of the model for the output file.
        distances (list[float]): Distances to calculate potential energy at.
        results (dict[str, dict[str, list[float | list[list[float]]]]]): Results dict
            to collect energies and forces at given distances for all diatomic curves.
            Will be updated in-place.

    Returns:
        dict[str, dict[str, list[float | list[list[float]]]]]: Potential energy and
            forces data for diatomic molecules.
    """
    # saving results in dict: {"symbol-symbol": {"energies": [...], "forces": [...]}}
    for idx, (z1, z2) in (pbar := tqdm(enumerate(pairs, start=1))):
        # Convert atomic numbers to symbols if needed
        elem1 = atom_num_symbol_map.get(z1, z1)
        elem2 = atom_num_symbol_map.get(z2, z2)
        formula = f"{elem1}-{elem2}"
        ef_dict = results.setdefault(formula, {"energies": [], "forces": []})

        len_e, len_f = len(ef_dict.get("energies", [])), len(ef_dict.get("forces", []))
        # skip if we have results for this formula and they match expected length
        if len_e == len_f == len(distances):
            continue

        pbar.set_description(
            f"{idx}/{len(pairs)} {formula} diatomic curve with {model_name}"
        )

        # reset ef_dict in case we had prior results
        results[formula] |= {"energies": [], "forces": []}
        try:
            for atoms in generate_diatomics(elem1, elem2, distances):
                results[formula]["energies"] += [calculator.get_potential_energy(atoms)]
                results[formula]["forces"] += [calculator.get_forces(atoms).tolist()]
        except Exception as exc:
            print(f"{idx}/{len(pairs)} {formula} failed: {exc}")
            continue

    return results
