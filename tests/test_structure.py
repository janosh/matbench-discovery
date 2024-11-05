import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.structure import (
    analyze_symmetry,
    perturb_structure,
    pred_vs_ref_struct_symmetry,
)


@pytest.mark.parametrize("gamma", [0.5, 1.0, 1.5, 2.0])
def test_perturb_structure(dummy_struct: Structure, gamma: float) -> None:
    perturbed = perturb_structure(dummy_struct, gamma=gamma)
    assert len(perturbed) == len(dummy_struct)

    for site, new in zip(dummy_struct, perturbed, strict=True):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # Check that different calls produce different results
    assert perturb_structure(dummy_struct, gamma=gamma) != perturb_structure(
        dummy_struct, gamma=gamma
    )


@pytest.fixture
def cubic_structure() -> Structure:
    lattice = Lattice.cubic(4.2)
    return Structure(lattice, ["Si", "Si"], [[0, 0, 0], [0.5, 0.5, 0.5]])


@pytest.fixture
def tetragonal_structure() -> Structure:
    lattice = Lattice.tetragonal(a=4.0, c=6.0)
    return Structure(
        lattice, ["Ti", "O", "O"], [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]]
    )


@pytest.fixture
def monoclinic_structure() -> Structure:
    lattice = Lattice.monoclinic(a=5.0, b=4.0, c=6.0, beta=100)
    return Structure(
        lattice,
        ["P", "O", "O", "O", "O"],
        [
            [0, 0, 0],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
        ],
    )


@pytest.mark.parametrize(
    "test_structure, expected_spg_num, expected_n_sym_ops",
    [
        ("cubic_structure", 229, 96),
        ("tetragonal_structure", 47, 8),
        ("monoclinic_structure", 3, 2),
    ],
)
def test_analyze_symmetry(
    request: pytest.FixtureRequest,
    test_structure: str,
    expected_spg_num: int,
    expected_n_sym_ops: int,
) -> None:
    structure = request.getfixturevalue(test_structure)
    df_sym = analyze_symmetry({"test_structure": structure})

    assert len(df_sym) == 1
    assert df_sym[Key.spg_num].iloc[0] == expected_spg_num
    assert df_sym[Key.n_sym_ops].iloc[0] == expected_n_sym_ops


def test_analyze_symmetry_multiple_structures(
    cubic_structure: Structure,
    tetragonal_structure: Structure,
    monoclinic_structure: Structure,
) -> None:
    structures = {
        "cubic": cubic_structure,
        "tetragonal": tetragonal_structure,
        "monoclinic": monoclinic_structure,
    }
    df_sym = analyze_symmetry(structures)

    assert len(df_sym) == 3
    assert list(df_sym[Key.spg_num]) == [229, 47, 3]
    assert list(df_sym[Key.n_sym_ops]) == [96, 8, 2]


def test_pred_vs_ref_struct_symmetry(
    cubic_structure: Structure, tetragonal_structure: Structure
) -> None:
    key = "structure"
    df_ml = analyze_symmetry({key: cubic_structure})
    df_dft = analyze_symmetry({key: tetragonal_structure})

    df_compared = pred_vs_ref_struct_symmetry(
        df_ml, df_dft, {key: cubic_structure}, {key: tetragonal_structure}
    )

    assert MbdKey.spg_num_diff in df_compared.columns
    assert MbdKey.n_sym_ops_diff in df_compared.columns
    assert df_compared[MbdKey.spg_num_diff].iloc[0] == 229 - 47
    n_sym_ops_ml, n_sym_ops_dft = 96, 8
    assert df_compared[MbdKey.n_sym_ops_diff].iloc[0] == n_sym_ops_ml - n_sym_ops_dft


def test_analyze_symmetry_perturbed_structure(cubic_structure: Structure) -> None:
    perturbed_structure = perturb_structure(cubic_structure, gamma=1.5)

    df_sym_original = analyze_symmetry({"original": cubic_structure})
    df_sym_perturbed = analyze_symmetry({"perturbed": perturbed_structure})

    assert not df_sym_original.equals(df_sym_perturbed)
    assert df_sym_perturbed[Key.spg_num].iloc[0] < df_sym_original[Key.spg_num].iloc[0]
    assert (
        df_sym_perturbed[Key.n_sym_ops].iloc[0] < df_sym_original[Key.n_sym_ops].iloc[0]
    )


def test_analyze_symmetry_supercell(cubic_structure: Structure) -> None:
    supercell = cubic_structure * (2, 2, 2)

    df_sym_original = analyze_symmetry({"original": cubic_structure})
    df_sym_supercell = analyze_symmetry({"supercell": supercell})

    assert list(df_sym_original) == list(df_sym_supercell)
    assert df_sym_original[Key.spg_num].iloc[0] == df_sym_supercell[Key.spg_num].iloc[0]


def test_pred_vs_ref_struct_symmetry_with_structures(
    cubic_structure: Structure, tetragonal_structure: Structure
) -> None:
    df_ml = pd.DataFrame({Key.structure: [cubic_structure]})
    df_dft = pd.DataFrame({Key.structure: [tetragonal_structure]})

    key = "struct1"
    df_ml_sym = analyze_symmetry({key: cubic_structure})
    df_dft_sym = analyze_symmetry({key: tetragonal_structure})

    df_ml = pd.concat([df_ml, df_ml_sym], axis=1)
    df_dft = pd.concat([df_dft, df_dft_sym], axis=1)

    df_ml_sym[Key.spg_num] - df_dft_sym[Key.spg_num]

    # must use same keys for both structures to match them in RMSD calculation
    df_compared = pred_vs_ref_struct_symmetry(
        df_ml, df_dft, {key: cubic_structure}, {key: tetragonal_structure}
    )

    assert set(df_compared) == {
        Key.choice_symbol,
        Key.hall_num,
        Key.hall_symbol,
        MbdKey.international_spg_name,
        Key.max_pair_dist,
        Key.n_rot_syms,
        Key.n_sym_ops,
        MbdKey.n_sym_ops_diff,
        Key.n_trans_syms,
        Key.point_group,
        Key.spg_num,
        MbdKey.spg_num_diff,
        Key.structure,
        MbdKey.structure_rmsd_vs_dft,
        Key.wyckoff_symbols,
    }


def test_analyze_symmetry_primitive_vs_conventional(cubic_structure: Structure) -> None:
    spg_analyzer = SpacegroupAnalyzer(cubic_structure)
    primitive_structure = spg_analyzer.get_primitive_standard_structure()

    key_conv, key_prim = "conventional", "primitive"
    df_conventional = analyze_symmetry({key_conv: cubic_structure})
    df_primitive = analyze_symmetry({key_prim: primitive_structure})
    assert df_conventional.index.name == Key.mat_id
    assert df_primitive.index.name == Key.mat_id
    assert df_primitive.index[0] == key_prim
    assert df_conventional.index[0] == key_conv

    cols_to_drop = [  # some columns differ between conventional and primitive structure
        Key.wyckoff_symbols,
        Key.n_sym_ops,
        Key.n_rot_syms,
        Key.n_trans_syms,
    ]
    df_primitive.index = df_conventional.index

    pd.testing.assert_frame_equal(
        df_conventional.drop(columns=cols_to_drop),
        df_primitive.drop(columns=cols_to_drop),
        check_index_type=False,
    )
