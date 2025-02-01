import pandas as pd
import pytest
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.structure import perturb_structure, symmetry


@pytest.mark.parametrize(
    "test_structure, expected_spg_num, expected_n_sym_ops",
    [
        ("cubic_struct", 229, 96),
        ("tetragonal_struct", 47, 8),
        ("monoclinic_struct", 3, 2),
    ],
)
def test_analyze_symmetry(
    request: pytest.FixtureRequest,
    test_structure: str,
    expected_spg_num: int,
    expected_n_sym_ops: int,
) -> None:
    structure = request.getfixturevalue(test_structure)
    df_sym = symmetry.get_sym_info_from_structs({"test_structure": structure})

    assert len(df_sym) == 1
    assert df_sym[Key.spg_num].iloc[0] == expected_spg_num
    assert df_sym[Key.n_sym_ops].iloc[0] == expected_n_sym_ops


def test_analyze_symmetry_multiple_structures(
    cubic_struct: Structure,
    tetragonal_struct: Structure,
    monoclinic_struct: Structure,
) -> None:
    structures = {
        "cubic": cubic_struct,
        "tetragonal": tetragonal_struct,
        "monoclinic": monoclinic_struct,
    }
    df_sym = symmetry.get_sym_info_from_structs(structures)

    assert len(df_sym) == 3
    assert list(df_sym[Key.spg_num]) == [229, 47, 3]
    assert list(df_sym[Key.n_sym_ops]) == [96, 8, 2]


def test_pred_vs_ref_struct_symmetry(
    cubic_struct: Structure, tetragonal_struct: Structure
) -> None:
    key = "structure"
    df_ml = symmetry.get_sym_info_from_structs({key: cubic_struct})
    df_dft = symmetry.get_sym_info_from_structs({key: tetragonal_struct})

    df_compared = symmetry.pred_vs_ref_struct_symmetry(
        df_ml, df_dft, {key: cubic_struct}, {key: tetragonal_struct}
    )

    assert MbdKey.spg_num_diff in df_compared.columns
    assert MbdKey.n_sym_ops_diff in df_compared.columns
    assert df_compared[MbdKey.spg_num_diff].iloc[0] == 229 - 47
    n_sym_ops_ml, n_sym_ops_dft = 96, 8
    assert df_compared[MbdKey.n_sym_ops_diff].iloc[0] == n_sym_ops_ml - n_sym_ops_dft


def test_analyze_symmetry_perturbed_structure(cubic_struct: Structure) -> None:
    perturbed_structure = perturb_structure(cubic_struct, gamma=1.5)

    df_sym_original = symmetry.get_sym_info_from_structs({"original": cubic_struct})
    df_sym_perturbed = symmetry.get_sym_info_from_structs(
        {"perturbed": perturbed_structure}
    )

    assert not df_sym_original.equals(df_sym_perturbed)
    assert df_sym_perturbed[Key.spg_num].iloc[0] < df_sym_original[Key.spg_num].iloc[0]
    assert (
        df_sym_perturbed[Key.n_sym_ops].iloc[0] < df_sym_original[Key.n_sym_ops].iloc[0]
    )


def test_analyze_symmetry_supercell(cubic_struct: Structure) -> None:
    supercell = cubic_struct * (2, 2, 2)

    df_sym_original = symmetry.get_sym_info_from_structs({"original": cubic_struct})
    df_sym_supercell = symmetry.get_sym_info_from_structs({"supercell": supercell})

    assert list(df_sym_original) == list(df_sym_supercell)
    assert df_sym_original[Key.spg_num].iloc[0] == df_sym_supercell[Key.spg_num].iloc[0]


def test_pred_vs_ref_struct_symmetry_with_structures(
    cubic_struct: Structure, tetragonal_struct: Structure
) -> None:
    df_ml = pd.DataFrame({Key.structure: [cubic_struct]})
    df_dft = pd.DataFrame({Key.structure: [tetragonal_struct]})

    key = "struct1"
    df_ml_sym = symmetry.get_sym_info_from_structs({key: cubic_struct})
    df_dft_sym = symmetry.get_sym_info_from_structs({key: tetragonal_struct})

    df_ml = pd.concat([df_ml, df_ml_sym], axis=1)
    df_dft = pd.concat([df_dft, df_dft_sym], axis=1)

    df_ml_sym[Key.spg_num] - df_dft_sym[Key.spg_num]

    # must use same keys for both structures to match them in RMSD calculation
    df_compared = symmetry.pred_vs_ref_struct_symmetry(
        df_ml, df_dft, {key: cubic_struct}, {key: tetragonal_struct}
    )

    assert set(df_compared) == {
        Key.hall_num,
        Key.hall_symbol,
        MbdKey.international_spg_name,
        Key.max_pair_dist,
        Key.n_rot_syms,
        Key.n_sym_ops,
        MbdKey.n_sym_ops_diff,
        Key.n_trans_syms,
        Key.spg_num,
        MbdKey.spg_num_diff,
        Key.structure,
        MbdKey.structure_rmsd_vs_dft,
        Key.wyckoff_symbols,
        Key.symprec,
        Key.angle_tolerance,
    }


def test_analyze_symmetry_primitive_vs_conventional(cubic_struct: Structure) -> None:
    spg_analyzer = SpacegroupAnalyzer(cubic_struct)
    primitive_structure = spg_analyzer.get_primitive_standard_structure()

    conventional_key, primitive_key = "conventional", "primitive"
    int_spg_key = "international_spg_name"
    df_conventional = symmetry.get_sym_info_from_structs(
        {conventional_key: cubic_struct}
    )
    df_primitive = symmetry.get_sym_info_from_structs(
        {primitive_key: primitive_structure}
    )
    assert df_conventional.index.name == Key.mat_id
    assert df_primitive.index.name == Key.mat_id
    assert df_primitive.index[0] == primitive_key
    assert df_conventional.index[0] == conventional_key
    assert df_conventional.index[0] == conventional_key
    assert df_conventional[int_spg_key].iloc[0] == ["m-3m", "m-3m"]
    assert df_primitive[int_spg_key].iloc[0] == ["m-3m"]

    cols_to_drop = [  # some columns differ between conventional and primitive structure
        Key.wyckoff_symbols,
        Key.n_sym_ops,
        Key.n_rot_syms,
        Key.n_trans_syms,
        int_spg_key,
    ]
    df_primitive.index = df_conventional.index

    pd.testing.assert_frame_equal(
        df_conventional.drop(columns=cols_to_drop),
        df_primitive.drop(columns=cols_to_drop),
        check_index_type=False,
    )


def test_analyze_symmetry_with_ase_atoms(cubic_struct: Structure) -> None:
    """Test analyze_symmetry with ASE Atoms objects."""
    from ase import Atoms
    from ase.spacegroup import crystal

    # Create an ASE Atoms object (simple cubic Cu)
    a = 3.6
    atoms = crystal("Cu", [(0, 0, 0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
    assert isinstance(atoms, Atoms)

    # Test with both Structure and Atoms
    df_ase = symmetry.get_sym_info_from_structs({"ase": atoms})

    assert len(df_ase) == 1
    assert df_ase.index.name == Key.mat_id
    assert isinstance(df_ase[Key.spg_num].iloc[0], int)
    assert isinstance(df_ase[Key.n_sym_ops].iloc[0], int)

    # Test mixed dictionary of Structure and Atoms
    df_mixed = symmetry.get_sym_info_from_structs({"pmg": cubic_struct, "ase": atoms})
    assert len(df_mixed) == 2
    assert df_mixed.index.tolist() == ["pmg", "ase"]
    assert all(isinstance(x, int) for x in df_mixed[Key.spg_num])
    assert all(isinstance(x, int) for x in df_mixed[Key.n_sym_ops])

    # should give same results for Structure and equivalent Atoms
    cubic_atoms = cubic_struct.to_ase_atoms()

    df_struct = symmetry.get_sym_info_from_structs({"struct": cubic_struct})
    df_atoms = symmetry.get_sym_info_from_structs({"atoms": cubic_atoms})

    # Reset index to avoid index name mismatch in frame comparison
    df_struct.index = df_atoms.index

    pd.testing.assert_frame_equal(df_struct, df_atoms)
