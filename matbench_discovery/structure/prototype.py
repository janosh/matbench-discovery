"""Get AFLOW prototype labels for crystal structures using Moyopy for symmetry
detection.
"""

import gzip
import itertools
import math
import os
import re
import string
from typing import Final

import ase
import yaml
from pymatgen.core import Composition
from pymatviz.typing import AnyStructure

module_dir = os.path.dirname(__file__)

with gzip.open(f"{module_dir}/wyckoff-position-multiplicities.yaml.gz") as gz_file:
    wyckoff_multiplicity_dict = yaml.safe_load(gz_file)

with gzip.open(f"{module_dir}/wyckoff-position-relabelings.yaml.gz") as gz_file:
    wyckoff_relabelings = yaml.safe_load(gz_file)

crys_sys_letters: Final[dict[str, str]] = {
    "Triclinic": "a",
    "Monoclinic": "m",
    "Orthorhombic": "o",
    "Tetragonal": "t",
    "Trigonal": "h",
    "Hexagonal": "h",
    "Cubic": "c",
}


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
        gcd = math.gcd(*map(int, composition.values()))
        reduced /= gcd

    return "".join(
        f"{elem}{int(amt) if abs(amt % 1) < amt_tol else amt}" if amt != 1 else elem
        for elem, amt in zip(
            string.ascii_uppercase,
            [reduced[str(key)] for key in sorted(reduced, key=str)],
            strict=False,
        )
    )


def canonicalize_wyckoffs(element_wyckoffs: str, spg_num: int) -> str:
    """Given an element ordering, canonicalize the associated Wyckoff positions
    based on the alphabetical weight of equivalent choices of origin.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores.
        spg_num (int | str): International space group number.

    Returns:
        str: element_wyckoff string with canonical ordering of the wyckoff letters.
    """
    scored_wyckoffs: list[tuple[str, int]] = []

    for trans in wyckoff_relabelings[spg_num]:
        wyckoffs = element_wyckoffs.translate(str.maketrans(trans))
        score: int = 0
        sorted_wyckoffs: list[str] = []

        for el_wyckoff_letter in wyckoffs.split("_"):
            # Split into alpha and numeric groups
            groups = [
                "".join(grp)
                for _, grp in itertools.groupby(el_wyckoff_letter, str.isalpha)
            ]
            letters = [char for char in groups if char.isalpha()]
            counts = [char for char in groups if char.isnumeric()]

            # Sort by Wyckoff letter and build string
            sorted_pairs = sorted(
                zip(counts, letters, strict=False), key=lambda x: x[1]
            )
            sorted_wyckoffs.append(
                "".join(
                    f"{count}{letter}" if count != "1" else letter
                    for count, letter in sorted_pairs
                )
            )
            score += sum(0 if letter == "A" else ord(letter) - 96 for letter in letters)

        scored_wyckoffs += [("_".join(sorted_wyckoffs), score)]

    return min(scored_wyckoffs, key=lambda x: (x[1], x[0]))[0]


def get_protostructure_label(
    struct: AnyStructure, *, symprec: float = 0.1, raise_errors: bool = False
) -> str | None:
    """Get AFLOW-style proto-structure label using Moyopy for symmetry detection.

    Args:
        struct (AnyStructure): pymatgen Structure object or ase.Atoms object.
        symprec (float): Initial symmetry precision for Moyopy. Defaults to 0.1.
        raise_errors (bool): Whether to raise ValueError for failing structures or
            return the error message as string instead of the prototype label. Defaults
            to False.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    import moyopy.interface
    from pymatgen.core import Structure

    pmg_struct: Structure
    if isinstance(struct, ase.Atoms):
        from pymatgen.io.ase import AseAtomsAdaptor

        converted = AseAtomsAdaptor.get_structure(struct)
        if not isinstance(converted, Structure):
            raise TypeError(f"Expected Structure, got {type(converted)}")
        pmg_struct = converted
    elif isinstance(struct, Structure):
        pmg_struct = struct
    else:
        raise TypeError(f"Expected Structure or Atoms, got {type(struct)}")

    moyo_cell = moyopy.interface.MoyoAdapter.from_py_obj(pmg_struct)
    symmetry_data = moyopy.MoyoDataset(moyo_cell, symprec=symprec)
    spg_num = symmetry_data.number

    # Group sites by orbit
    orbit_groups: dict[int, list[int]] = {}

    for idx, orbit_id in enumerate(symmetry_data.orbits):
        if orbit_id not in orbit_groups:
            orbit_groups[orbit_id] = []
        orbit_groups[orbit_id].append(idx)

    # Create equivalent_wyckoff_labels from orbit groups
    element_dict: dict[str, int] = {}
    element_wyckoffs: list[str] = []

    equivalent_wyckoff_labels: list[tuple[str, int, str]] = [
        (
            pmg_struct.species[orbit[0]].symbol,
            len(orbit),
            symmetry_data.wyckoffs[orbit[0]].translate(
                str.maketrans("", "", string.digits)
            ),
        )
        for orbit in orbit_groups.values()
    ]
    equivalent_wyckoff_labels = sorted(
        equivalent_wyckoff_labels, key=lambda x: (x[0], x[2])
    )

    for el, group in itertools.groupby(equivalent_wyckoff_labels, key=lambda x: x[0]):
        # NOTE create a list from the iterator so that we can use it without exhausting
        label_group = list(group)
        element_dict[el] = sum(
            wyckoff_multiplicity_dict[spg_num][itm[2]] for itm in label_group
        )
        wyckoff_counts_str = "".join(
            f"{len(list(occurrences))}{wyk_letter}"
            for wyk_letter, occurrences in itertools.groupby(
                label_group, key=lambda x: x[2]
            )
        )

        element_wyckoffs += [wyckoff_counts_str]

    all_wyckoffs = canonicalize_wyckoffs("_".join(element_wyckoffs), spg_num)

    # Build prototype label
    prototype_label = (
        f"{get_prototype_formula(pmg_struct.composition)}_"
        f"{symmetry_data.pearson_symbol}_"
        f"{spg_num}_{all_wyckoffs}:{pmg_struct.chemical_system}"
    )

    # Verify multiplicities match composition
    observed_formula = Composition(element_dict).reduced_formula
    if observed_formula != pmg_struct.reduced_formula:
        err_msg = (
            f"Invalid WP multiplicities - {prototype_label}, expected "
            f"{observed_formula} to be {pmg_struct.reduced_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)
        return err_msg

    return prototype_label


def count_wyckoff_positions(protostructure_label: str) -> int:
    """Count number of Wyckoff positions in an Aflow-style protostructure label.

    Adapted from https://github.com/CompRhys/aviary/blob/181e2b2b2/aviary/wren/utils.py#L472.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: number of distinct Wyckoff positions in protostructure_label

    Raises:
        ValueError: On invalid protostructure labels
    """
    aflow_label, *_chem_sys = protostructure_label.split(":")  # remove chemical system
    # discard prototype formula and spg symbol and spg number
    parts = aflow_label.split("_", maxsplit=3)
    if len(parts) < 4:
        raise ValueError("Invalid protostructure label: missing required parts")

    # throw Wyckoff positions for all elements together
    wyk_letters = parts[-1].replace("_", "")
    if not wyk_letters.strip():
        raise ValueError("Invalid protostructure label: empty Wyckoff positions")

    wyk_list = re.split("[A-z]", wyk_letters)[:-1]  # split on every letter

    # count 1 for letters without prefix
    return sum(1 if len(x) == 0 else int(x) for x in wyk_list)
