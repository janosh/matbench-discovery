"""Get AFLOW prototype labels for crystal structures using Moyopy for symmetry
detection.
"""

import gzip
import itertools
import math
import os
import re
import string

import ase
import yaml
from pymatgen.core import Composition
from pymatviz.typing import AnyStructure

module_dir = os.path.dirname(__file__)

with gzip.open(f"{module_dir}/wyckoff-position-multiplicities.yaml.gz") as gz_file:
    wyckoff_multiplicity_dict = yaml.safe_load(gz_file)

with gzip.open(f"{module_dir}/wyckoff-position-relabelings.yaml.gz") as gz_file:
    wyckoff_relabelings = yaml.safe_load(gz_file)


def get_prototype_formula(composition: Composition, amt_tol: float = 1e-8) -> str:
    """Get anonymized formula for a Composition where species are in alphabetical
    order and assigned ascending letters. This format is used in the Aflow
    structure prototype labeling scheme.

    Args:
        composition (Composition): Pymatgen Composition to process
        amt_tol (float): Tolerance for amount of species. Defaults to 1e-8.

    Returns:
        str: anonymized formula where the species are in alphabetical order
    """
    reduced = composition.element_composition
    if all(elem_amount == int(elem_amount) for elem_amount in composition.values()):
        reduced /= math.gcd(*map(int, composition.values()))

    return "".join(
        f"{elem}{int(amt) if abs(amt % 1) < amt_tol else amt}" if amt != 1 else elem
        for elem, amt in zip(
            string.ascii_uppercase,
            [reduced[key] for key in sorted(reduced, key=str)],
            strict=False,
        )
    )


def canonicalize_wyckoffs(element_wyckoffs: str, spg_num: int) -> str:
    """Given an element ordering, canonicalize the associated Wyckoff positions
    based on the alphabetical weight of equivalent choices of origin.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores, each
            letter preceded by its count (e.g. "1a2c_4b").
        spg_num (int | str): International space group number.

    Returns:
        str: element_wyckoff string with canonical ordering of the wyckoff letters.
    """
    scored_wyckoffs: list[tuple[int, str]] = []  # (score, relabeled wyckoffs)

    for trans in wyckoff_relabelings[spg_num]:
        score = 0
        sorted_wyckoffs: list[str] = []
        for el_wyckoffs in element_wyckoffs.translate(trans).split("_"):
            count_letter_pairs = re.findall(r"(\d+)([a-zA-Z])", el_wyckoffs)
            sorted_wyckoffs.append(
                "".join(
                    f"{count}{letter}" if count != "1" else letter
                    for count, letter in sorted(
                        count_letter_pairs, key=lambda pair: pair[1]
                    )
                )
            )
            score += sum(
                0 if letter == "A" else ord(letter) - 96
                for _, letter in count_letter_pairs
            )
        scored_wyckoffs.append((score, "_".join(sorted_wyckoffs)))

    return min(scored_wyckoffs)[1]


def get_protostructure_label(struct: AnyStructure, *, symprec: float = 0.1) -> str:
    """Get AFLOW-style proto-structure label using Moyopy for symmetry detection.

    Args:
        struct (AnyStructure): pymatgen Structure object or ase.Atoms object.
        symprec (float): Initial symmetry precision for Moyopy. Defaults to 0.1.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys`.

    Raises:
        ValueError: If the detected Wyckoff multiplicities do not reproduce the
            structure's composition (failed symmetry detection).
    """
    import moyopy.interface
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

    pmg_struct = (
        AseAtomsAdaptor.get_structure(struct)
        if isinstance(struct, ase.Atoms)
        else struct
    )
    if not isinstance(pmg_struct, Structure):
        raise TypeError(f"Expected Structure or Atoms, got {type(struct)}")

    moyo_cell = moyopy.interface.MoyoAdapter.from_py_obj(pmg_struct)
    symmetry_data = moyopy.MoyoDataset(moyo_cell, symprec=symprec)
    spg_num = symmetry_data.number

    # one (element, Wyckoff letter) per orbit, sorted: all its sites share both
    site_per_orbit = {orbit: idx for idx, orbit in enumerate(symmetry_data.orbits)}
    orbit_labels = sorted(
        (
            pmg_struct.species[idx].symbol,
            symmetry_data.wyckoffs[idx].lstrip(string.digits),
        )
        for idx in site_per_orbit.values()
    )
    element_dict: dict[str, int] = {}
    element_wyckoffs: list[str] = []
    for element, group in itertools.groupby(orbit_labels, key=lambda label: label[0]):
        letters = [letter for _, letter in group]
        element_dict[element] = sum(
            wyckoff_multiplicity_dict[spg_num][letter] for letter in letters
        )
        element_wyckoffs.append(
            "".join(
                f"{len(list(same))}{letter}"
                for letter, same in itertools.groupby(letters)
            )
        )

    all_wyckoffs = canonicalize_wyckoffs("_".join(element_wyckoffs), spg_num)

    prototype_label = (
        f"{get_prototype_formula(pmg_struct.composition)}_"
        f"{symmetry_data.pearson_symbol}_"
        f"{spg_num}_{all_wyckoffs}:{pmg_struct.chemical_system}"
    )

    # Verify multiplicities match composition
    observed_formula = Composition(element_dict).reduced_formula
    if observed_formula != pmg_struct.reduced_formula:
        raise ValueError(
            f"Invalid WP multiplicities - {prototype_label}, expected "
            f"{observed_formula} to be {pmg_struct.reduced_formula}"
        )

    return prototype_label
