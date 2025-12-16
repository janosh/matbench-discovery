from typing import Final

import pytest
from pymatgen.core.structure import Lattice, Structure

from matbench_discovery import TEST_FILES
from matbench_discovery.structure import prototype

NaCl = Structure(
    lattice=[[2, 2, 0], [0, 2, 2], [2, 0, 2]],
    species=["Na", "Cl"],
    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
)
CsCl = Structure(
    lattice=[[4, 0, 0], [0, 4, 0], [0, 0, 4]],
    species=["Cs", "Cl"],
    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
)
ZnO_zincblende = Structure(
    lattice=[[2, 2, 0], [0, 2, 2], [2, 0, 2]],
    species=["Zn", "O"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)
ZnO_wurtzite = Structure(
    lattice=Lattice.from_parameters(
        a=3.8227, b=3.8227, c=6.2607, alpha=90, beta=90, gamma=120
    ),
    species=["Zn", "O", "Zn", "O"],
    coords=[
        [1 / 3, 2 / 3, 0],
        [2 / 3, 1 / 3, 0.3748],
        [2 / 3, 1 / 3, 1 / 2],
        [1 / 3, 2 / 3, 1 / 2 + 0.3748],
    ],
)
CaF2_fluorite = Structure(
    lattice=Lattice.from_parameters(a=3.9, b=3.9, c=3.9, alpha=60, beta=60, gamma=60),
    species=["Ca", "F", "F"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
)
Cu3Au = Structure(
    lattice=[[3.75, 0, 0], [0, 3.75, 0], [0, 0, 3.75]],
    species=["Au", "Cu", "Cu", "Cu"],
    coords=[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
)
Fe3Al_DO3 = Structure(
    lattice=[[5.76, 0, 0], [0, 5.76, 0], [0, 0, 5.76]],
    species=["Al", "Fe", "Fe", "Fe", "Al", "Fe", "Fe", "Fe"],
    coords=[
        [0, 0, 0],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0],
        [0.75, 0.75, 0.25],
        [0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0, 0.5],
        [0.75, 0.25, 0.75],
    ],
)
SrTiO3_perovskite = Structure(
    lattice=[[3.9, 0, 0], [0, 3.9, 0], [0, 0, 3.9]],
    species=["Sr", "Ti", "O", "O", "O"],
    coords=[[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
)
esseneite = Structure.from_file(f"{TEST_FILES}/structures/esseneite.cif.gz")

# Labels from moyopy==0.3.1 and spglib==2.5.0 differed for a few cases which were
# tracked and one by one resolved in https://github.com/CompRhys/aviary/pull/96
C20Cd8H14N4O4S_oP800 = Structure.from_file(
    f"{TEST_FILES}/structures/A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c@C-Cd-H-N-O-S.POSCAR.xz"
)
BaTiO3_perovskite = Structure.from_file(
    f"{TEST_FILES}/structures/AB3C_cP5_221_a_c_b@Ba-O-Ti.POSCAR.xz"
)

TEST_CASES: Final[tuple[tuple[Structure, str], ...]] = (
    (NaCl, "AB_cF8_225_a_b:Cl-Na"),
    (CsCl, "AB_cP2_221_a_b:Cl-Cs"),
    (ZnO_zincblende, "AB_cF8_216_a_c:O-Zn"),
    (ZnO_wurtzite, "AB_hP4_186_b_b:O-Zn"),
    (CaF2_fluorite, "AB2_cF12_225_a_c:Ca-F"),
    (Cu3Au, "AB3_cP4_221_a_c:Au-Cu"),
    (Fe3Al_DO3, "AB3_tP4_115_a_cg:Al-Fe"),
    (SrTiO3_perovskite, "A3BC_cP5_221_c_a_b:O-Sr-Ti"),
    (esseneite, "ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si"),
    (C20Cd8H14N4O4S_oP800, "A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c:C-Cd-H-N-O-S"),
    (BaTiO3_perovskite, "AB3C_cP5_221_a_c_b:Ba-O-Ti"),
)


@pytest.mark.parametrize("struct, label", TEST_CASES)
def test_get_protostructure_label(struct: Structure, label: str) -> None:
    """Check that moyopy gives correct protostructure label for simple cases."""
    moyopy_label = prototype.get_protostructure_label(struct)
    assert moyopy_label == label, f"{moyopy_label=} != {label}"


def test_get_protostructure_label_with_ase_atoms() -> None:
    """Test that get_protostructure_label() works with ASE Atoms input."""
    from ase import Atoms
    from pymatgen.io.ase import AseAtomsAdaptor

    # Convert NaCl structure to ASE Atoms
    NaCl_atoms = AseAtomsAdaptor.get_atoms(NaCl)
    assert isinstance(NaCl_atoms, Atoms)
    # Get prototype label for both Structure and Atoms versions
    struct_label = prototype.get_protostructure_label(NaCl)
    atoms_label = prototype.get_protostructure_label(NaCl_atoms)

    # Check that both give the same result
    assert atoms_label == struct_label == "AB_cF8_225_a_b:Cl-Na"

    # Test with a more complex structure
    SrTiO3_atoms = AseAtomsAdaptor.get_atoms(SrTiO3_perovskite)
    assert isinstance(SrTiO3_atoms, Atoms)
    atoms_label = prototype.get_protostructure_label(SrTiO3_atoms)
    assert atoms_label == "A3BC_cP5_221_c_a_b:O-Sr-Ti"


@pytest.mark.parametrize(
    "protostructure_label, expected",
    [
        ("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si", 6),  # esseneite
        ("A6B11CD7_aP50_2_6i_ac10i_i_7i:C-H-N-O", 26),
        ("foo_bar_47_abc_A_b:X-Y-Z", 5),
        # Simple cases
        ("AB_cF8_225_a_b:Na-Cl", 2),  # NaCl structure
        ("AB_cP2_221_a_b:Cs-Cl", 2),  # CsCl structure
        ("AB2_cF12_225_a_c:Ca-F", 2),  # Fluorite structure
        # Multiple letters without numbers
        ("AB3_cP4_221_a_c:Au-Cu", 2),  # Cu3Au structure
        ("AB3_tP4_115_a_cg:Al-Fe", 3),  # DO3 structure
        # Complex cases with multiple numbers
        ("A2B3C4_oP18_62_2c_3d_4a:Fe-O-Si", 9),
        ("A3B2C6D_hP24_194_3a_2b_6c_d:Al-Fe-O-Si", 12),
        # Edge cases
        ("A_cF4_225_a:Cu", 1),  # Single Wyckoff position
        ("ABCD_cP4_221_abcd:W-X-Y-Z", 4),  # All single letters
        ("A10B10C10_cF120_225_10a_10b_10c:Fe-Ni-Al", 30),  # Large numbers
        ("A2B2_mP8_10_2a_2a:Si-O", 4),  # Repeated positions
    ],
)
def test_count_wyckoff_positions(protostructure_label: str, expected: int) -> None:
    """Test counting Wyckoff positions in prototype labels with various patterns."""
    count = prototype.count_wyckoff_positions(protostructure_label)
    assert isinstance(count, int)
    assert count == expected, (
        f"{expected=} Wyckoff positions for {protostructure_label=}, got {count=}"
    )


# Test error cases
@pytest.mark.parametrize(
    "invalid_label, err_msg",
    [
        # Empty Wyckoff positions before colon
        ("AB_cF8_225_:Na-Cl", "Invalid protostructure label: empty Wyckoff positions"),
        # Empty Wyckoff positions before colon
        ("AB_cF8_225:-Cl", "Invalid protostructure label: missing required parts"),
        # jumbled mess
        ("AB_cF8asx-225_a+:", "Invalid protostructure label: missing required parts"),
    ],
)
def test_count_wyckoff_positions_errors(invalid_label: str, err_msg: str) -> None:
    """Test error handling for invalid prototype labels.

    The function raises ValueError when:
    - The label doesn't have enough underscore-separated parts
    - It can't find Wyckoff positions
    - The Wyckoff position part is empty
    - The Wyckoff position part before colon is empty
    - The chemical system (after colon) is missing
    """
    with pytest.raises(ValueError, match=err_msg):
        prototype.count_wyckoff_positions(invalid_label)
