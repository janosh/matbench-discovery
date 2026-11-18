"""Regression test for the Aflow-style protostructure labels behind WBM's uniq-proto
filter, now computed by matscikit.protostructure_label.
"""

from collections.abc import Sequence

import numpy as np
import pytest
from matscikit import Structure, protostructure_label
from matscikit.interop import structure_from_pymatgen_dict, structure_to_pymatgen_dict
from matscikit.io import read_poscar

from matbench_discovery import TEST_FILES

FCC_PRIMITIVE = [[2, 2, 0], [0, 2, 2], [2, 0, 2]]


def lattice_from_parameters(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> np.ndarray:
    """Row lattice vectors from cell lengths (Å) and angles (degrees), a along x."""
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    vec_b = b * np.array([np.cos(gamma), np.sin(gamma), 0])
    c_x = c * np.cos(beta)
    c_y = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    return np.array([[a, 0, 0], vec_b, [c_x, c_y, np.sqrt(c**2 - c_x**2 - c_y**2)]])


def frac_structure(
    lattice: object, species: list[str], coords: Sequence[Sequence[float]]
) -> Structure:
    """Structure from fractional coordinates."""
    return Structure(lattice, species, coords, fractional=True)


def poscar_structure(name: str) -> Structure:
    """Structure of a POSCAR test file."""
    return read_poscar(f"{TEST_FILES}/structures/{name}.POSCAR.xz").structure


fcc_corner = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
TEST_CASES: dict[str, tuple[Structure, str]] = {
    "NaCl": (
        frac_structure(FCC_PRIMITIVE, ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
        "AB_cF8_225_a_b:Cl-Na",
    ),
    "CsCl": (
        frac_structure(4 * np.eye(3), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
        "AB_cP2_221_a_b:Cl-Cs",
    ),
    "ZnO-zincblende": (
        frac_structure(FCC_PRIMITIVE, ["Zn", "O"], [[0, 0, 0], [0.25, 0.25, 0.25]]),
        "AB_cF8_216_a_c:O-Zn",
    ),
    "ZnO-wurtzite": (
        frac_structure(
            lattice_from_parameters(3.8227, 3.8227, 6.2607, 90, 90, 120),
            ["Zn", "O", "Zn", "O"],
            [
                [1 / 3, 2 / 3, 0],
                [2 / 3, 1 / 3, 0.3748],
                [2 / 3, 1 / 3, 1 / 2],
                [1 / 3, 2 / 3, 1 / 2 + 0.3748],
            ],
        ),
        "AB_hP4_186_b_b:O-Zn",
    ),
    "CaF2-fluorite": (
        frac_structure(
            lattice_from_parameters(3.9, 3.9, 3.9, 60, 60, 60),
            ["Ca", "F", "F"],
            [[0, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
        ),
        "AB2_cF12_225_a_c:Ca-F",
    ),
    "Cu3Au": (
        frac_structure(3.75 * np.eye(3), ["Au", "Cu", "Cu", "Cu"], fcc_corner),
        "AB3_cP4_221_a_c:Au-Cu",
    ),
    "Fe3Al-DO3": (
        frac_structure(
            5.76 * np.eye(3),
            ["Al", "Fe", "Fe", "Fe"] * 2,
            [
                [0, 0, 0],
                [0.25, 0.25, 0.25],
                [0.5, 0.5, 0],
                [0.75, 0.75, 0.25],
                [0, 0.5, 0.5],
                [0.25, 0.75, 0.75],
                [0.5, 0, 0.5],
                [0.75, 0.25, 0.75],
            ],
        ),
        "AB3_tP4_115_a_cg:Al-Fe",
    ),
    "SrTiO3-perovskite": (
        frac_structure(
            3.9 * np.eye(3),
            ["Sr", "Ti", "O", "O", "O"],
            [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        ),
        "A3BC_cP5_221_c_a_b:O-Sr-Ti",
    ),
    "esseneite": (
        poscar_structure("esseneite"),
        "ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si",
    ),
    # labels from moyopy and spglib once differed for these two, resolved in
    # https://github.com/CompRhys/aviary/pull/96
    "C20Cd8H14N4O4S": (
        poscar_structure("A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c@C-Cd-H-N-O-S"),
        "A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c:C-Cd-H-N-O-S",
    ),
    "BaTiO3-perovskite": (
        poscar_structure("AB3C_cP5_221_a_c_b@Ba-O-Ti"),
        "AB3C_cP5_221_a_c_b:Ba-O-Ti",
    ),
}


@pytest.mark.parametrize("as_dict", [False, True])
@pytest.mark.parametrize(("struct", "label"), TEST_CASES.values(), ids=TEST_CASES)
def test_protostructure_label(struct: Structure, label: str, as_dict: bool) -> None:
    """Labels match the deleted matbench implementation, also via a pymatgen dict."""
    if as_dict:
        struct = structure_from_pymatgen_dict(structure_to_pymatgen_dict(struct))
    assert protostructure_label(struct) == label
