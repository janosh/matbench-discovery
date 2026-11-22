"""Tests for crystal symmetry analysis helpers."""

import numpy as np
import pandas as pd
import pytest
from kemio import Structure, symmetry_dataset
from kemio.interop import structure_to_pymatgen_dict
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.structure import symmetry


def test_analyze_symmetry_multiple_structures(
    cubic_struct: Structure,
    tetragonal_struct: Structure,
    monoclinic_struct: Structure,
) -> None:
    """Space group numbers, symbols and operation counts for three lattice types."""
    structures = {
        "cubic": cubic_struct,
        "tetragonal": tetragonal_struct,
        "monoclinic": monoclinic_struct,
    }
    df_sym = symmetry.get_sym_info_from_structs(structures)

    assert len(df_sym) == 3
    assert list(df_sym[Key.spg_num]) == [229, 47, 3]
    # output schema of the geo-opt symmetry CSVs
    assert list(df_sym) == [
        Key.spg_num,
        Key.hall_num,
        MbdKey.international_spg_name,
        Key.wyckoff_symbols,
        Key.n_sym_ops,
        Key.n_rot_syms,
        Key.n_trans_syms,
        Key.hall_symbol,
        Key.symprec,
        Key.angle_tolerance,
    ]
    for col in (Key.n_sym_ops, Key.n_rot_syms, Key.n_trans_syms):
        assert list(df_sym[col]) == [96, 8, 2]
    assert list(df_sym[MbdKey.international_spg_name]) == ["I m -3 m", "P m m m", "P 2"]
    # Hall symbols are the setting-specific symbols, not the international short names
    assert df_sym[Key.hall_symbol].iloc[0] == "-I 4 2 3"
    assert df_sym[Key.hall_num].iloc[0] == 529


def test_pred_vs_ref_struct_symmetry(tetragonal_struct: Structure) -> None:
    """Symmetry diffs use ML minus DFT, RMSD is measured, and columns survive."""
    key = "structure"
    # displace one O off its Wyckoff site: lowers symmetry and gives a nonzero RMSD
    positions = tetragonal_struct.positions.copy()
    positions[1] += [0.05, 0.02, 0.01]
    pred_struct = tetragonal_struct.replace(positions=positions)
    df_ml = symmetry.get_sym_info_from_structs({key: pred_struct})
    df_dft = symmetry.get_sym_info_from_structs({key: tetragonal_struct})
    df_ml_in = df_ml.assign(**{Key.structure: [pred_struct]})
    df_dft_in = df_dft.assign(**{Key.structure: [tetragonal_struct]})

    # must use same keys for both structures to match them in RMSD calculation
    cols_in = list(df_ml_in)
    df_compared = symmetry.pred_vs_ref_struct_symmetry(
        df_ml_in, df_dft_in, {key: pred_struct}, {key: tetragonal_struct}
    )
    assert list(df_ml_in) == cols_in  # returns a copy, input left unchanged

    spg_ml, spg_dft = df_ml[Key.spg_num].iloc[0], df_dft[Key.spg_num].iloc[0]
    assert spg_dft == 47
    assert spg_ml != spg_dft
    assert df_compared[MbdKey.spg_num_diff].iloc[0] == spg_ml - spg_dft
    n_ops_ml, n_ops_dft = df_ml[Key.n_sym_ops].iloc[0], df_dft[Key.n_sym_ops].iloc[0]
    assert n_ops_ml < n_ops_dft == 8
    assert df_compared[MbdKey.n_sym_ops_diff].iloc[0] == n_ops_ml - n_ops_dft
    rmsd = df_compared[MbdKey.structure_rmsd_vs_dft].iloc[0]
    assert rmsd is not None
    assert 0 < rmsd < 0.1

    # IDs without a DFT reference keep their row with missing RMSD and max distance
    preds = {key: pred_struct, "pred-only": tetragonal_struct}
    df_partial = symmetry.pred_vs_ref_struct_symmetry(
        symmetry.get_sym_info_from_structs(preds),
        df_dft,
        preds,
        {key: tetragonal_struct},
    )
    rms_cols = [MbdKey.structure_rmsd_vs_dft, Key.max_pair_dist]
    assert (df_partial[rms_cols].dtypes == "float64").all()
    assert df_partial.loc["pred-only", rms_cols].isna().all()
    assert (
        df_partial.loc[key, rms_cols].tolist()
        == df_compared.loc[key, rms_cols].tolist()
    )
    assert set(df_compared) == {
        Key.hall_num,
        Key.hall_symbol,
        MbdKey.international_spg_name,
        Key.max_pair_dist,
        Key.n_sym_ops,
        Key.n_rot_syms,
        Key.n_trans_syms,
        MbdKey.n_sym_ops_diff,
        Key.spg_num,
        MbdKey.spg_num_diff,
        Key.structure,
        MbdKey.structure_rmsd_vs_dft,
        Key.wyckoff_symbols,
        Key.symprec,
        Key.angle_tolerance,
    }


def test_geo_opt_relaxed_boundary_cell() -> None:
    """A relaxed WBM lattice on a Niggli boundary must not terminate the RMSD batch."""
    lattice = [
        [3.5162634007, 1.7628e-6, -1.02026e-5],
        [-5.28446e-5, 8.4037724935, 4.8518887469],
        [-2.1294e-6, 3.8997e-6, 9.7038417713],
    ]
    structures = {
        "boundary": Structure(
            lattice, ["Si", "O"], [[0, 0, 0], [0.25, 0.26, 0.2]], fractional=True
        )
    }
    summary = symmetry.get_sym_info_from_structs(structures)
    result = symmetry.pred_vs_ref_struct_symmetry(
        summary, summary, structures, structures
    )
    assert result.loc["boundary", MbdKey.structure_rmsd_vs_dft] == 0.0


@pytest.mark.parametrize("cell", ["supercell", "primitive"])
def test_analyze_symmetry_is_cell_choice_invariant(
    cubic_struct: Structure, cell: str
) -> None:
    """Supercells and primitive cells keep the conventional cell's space group."""
    if cell == "supercell":
        shifts = np.array(
            [[i, j, k] for i in range(2) for j in range(2) for k in range(2)]
        )
        frac = ((cubic_struct.frac_coords[None] + shifts[:, None]) / 2).reshape(-1, 3)
        other_struct = Structure(
            2 * cubic_struct.lattice, cubic_struct.symbols * 8, frac, fractional=True
        )
    else:
        sym = symmetry_dataset(cubic_struct)
        symbols = ["Si"] * len(sym.prim_std_numbers)
        other_struct = Structure(
            sym.prim_std_lattice,
            symbols,
            sym.prim_std_frac_coords,
            fractional=True,
        )
    df_conventional = symmetry.get_sym_info_from_structs({"conventional": cubic_struct})
    df_other = symmetry.get_sym_info_from_structs({cell: other_struct})
    for df_sym, mat_id in ((df_conventional, "conventional"), (df_other, cell)):
        assert df_sym.index.name == Key.mat_id
        assert df_sym.index.tolist() == [mat_id]
        # space group name is a property of the lattice, not of the cell choice
        assert df_sym[MbdKey.international_spg_name].iloc[0] == "I m -3 m"

    cols_to_drop = [  # columns that depend on the cell choice
        Key.wyckoff_symbols,
        Key.n_sym_ops,
        Key.n_rot_syms,
        Key.n_trans_syms,
    ]
    df_other.index = df_conventional.index
    pd.testing.assert_frame_equal(
        df_conventional.drop(columns=cols_to_drop),
        df_other.drop(columns=cols_to_drop),
        check_index_type=False,
    )


def test_analyze_symmetry_accepts_structure_dicts(cubic_struct: Structure) -> None:
    """Structure dicts give the right symmetry and match equivalent Structures."""
    fcc_coords = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    fcc_cu = Structure(3.6 * np.eye(3), ("Cu",) * 4, fcc_coords, fractional=True)
    fcc_dict = structure_to_pymatgen_dict(fcc_cu)

    df_fcc = symmetry.get_sym_info_from_structs({"fcc": fcc_dict})

    assert df_fcc.index.tolist() == ["fcc"]
    assert df_fcc.index.name == Key.mat_id
    fcc_row = df_fcc.loc["fcc"]
    # conventional 4-atom fcc cell: 48 point ops x 4 centering translations
    assert fcc_row[Key.spg_num] == 225
    assert fcc_row[MbdKey.international_spg_name] == "F m -3 m"
    assert fcc_row[Key.hall_symbol] == "-F 4 2 3"
    assert fcc_row[Key.n_sym_ops] == 192
    assert fcc_row[Key.wyckoff_symbols] == ["a"] * 4

    df_mixed = symmetry.get_sym_info_from_structs(
        {"obj": cubic_struct, "dict": fcc_dict}
    )
    assert df_mixed.index.tolist() == ["obj", "dict"]
    assert df_mixed[Key.spg_num].tolist() == [229, 225]

    # should give same results for Structure and equivalent dict
    df_struct = symmetry.get_sym_info_from_structs({"struct": cubic_struct})
    df_dict = symmetry.get_sym_info_from_structs(
        {"struct": structure_to_pymatgen_dict(cubic_struct)}
    )
    pd.testing.assert_frame_equal(df_struct, df_dict)
    with pytest.raises(TypeError, match=r"got builtins.int for .bad."):
        symmetry.get_sym_info_from_structs({"bad": 42})  # ty: ignore[invalid-argument-type]


def test_analyze_symmetry_failures(cubic_struct: Structure) -> None:
    """Failed symmetry searches give NaN rows; empty structures raise naming the ID."""
    overlapping = Structure(4 * np.eye(3), ["Na", "Cl"], np.zeros((2, 3)))
    df_sym = symmetry.get_sym_info_from_structs(
        {"ok": cubic_struct, "overlapping": overlapping}
    )
    assert df_sym[Key.spg_num].tolist()[0] == 229
    cols = [Key.spg_num, Key.n_sym_ops, Key.hall_symbol, Key.wyckoff_symbols]
    assert df_sym.loc["overlapping", cols].isna().all()
    # integer columns stay integer next to the failed row (no 229.0 in CSVs)
    int_cols = [Key.spg_num, Key.hall_num, Key.n_sym_ops, Key.n_rot_syms]
    int_cols += [Key.n_trans_syms]
    assert (df_sym[int_cols].dtypes == "Int64").all()
    csv_rows = df_sym[int_cols].to_csv(header=False).splitlines()
    assert csv_rows == ["ok,229,529,96,96,96", "overlapping,,,,,"]
    empty = Structure(4 * np.eye(3), [], np.zeros((0, 3)))
    with pytest.raises(ValueError, match="structure 1: ") as exc_info:
        symmetry.get_sym_info_from_structs({"ok": cubic_struct, "empty": empty})
    assert exc_info.value.__notes__ == [f"{Key.mat_id}='empty'"]
    with pytest.raises(ValueError, match="symprec"):
        symmetry.get_sym_info_from_structs({}, symprec=-1)
