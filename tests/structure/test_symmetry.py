"""Tests for crystal symmetry analysis helpers."""

import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from matbench_discovery.enums import Key
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
        Key.international_spg_name,
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
    assert list(df_sym[Key.international_spg_name]) == ["I m -3 m", "P m m m", "P 2"]
    # Hall symbols are the setting-specific symbols, not the international short names
    assert df_sym[Key.hall_symbol].iloc[0] == "-I 4 2 3"
    assert df_sym[Key.hall_num].iloc[0] == 529


def test_pred_vs_ref_struct_symmetry(tetragonal_struct: Structure) -> None:
    """Symmetry diffs use ML minus DFT, RMSD is measured, and columns survive."""
    key = "structure"
    pred_struct = tetragonal_struct.copy()
    # displace one O off its Wyckoff site: lowers symmetry and gives a nonzero RMSD
    pred_struct.translate_sites([1], [0.05, 0.02, 0.01], frac_coords=False)
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
    assert df_compared[Key.spg_num_diff].iloc[0] == spg_ml - spg_dft
    n_ops_ml, n_ops_dft = (
        df_ml[Key.n_sym_ops].iloc[0],
        df_dft[Key.n_sym_ops].iloc[0],
    )
    assert n_ops_ml < n_ops_dft == 8
    assert df_compared[Key.n_sym_ops_diff].iloc[0] == n_ops_ml - n_ops_dft
    assert 0 < df_compared[Key.structure_rmsd_vs_dft].iloc[0] < 0.1

    # IDs without a DFT reference keep their row with missing RMSD and max distance
    preds = {key: pred_struct, "pred-only": tetragonal_struct}
    df_partial = symmetry.pred_vs_ref_struct_symmetry(
        symmetry.get_sym_info_from_structs(preds),
        df_dft,
        preds,
        {key: tetragonal_struct},
    )
    rms_cols = [Key.structure_rmsd_vs_dft, Key.max_pair_dist]
    assert (df_partial[rms_cols].dtypes == "float64").all()
    assert df_partial.loc["pred-only", rms_cols].isna().all()
    assert (
        df_partial.loc[key, rms_cols].tolist()
        == df_compared.loc[key, rms_cols].tolist()
    )
    # all input columns survive in order, followed by the compared ones
    diff_cols = [Key.spg_num_diff, Key.n_sym_ops_diff]
    assert list(df_compared) == [*cols_in, *diff_cols, *rms_cols]


def test_geo_opt_relaxed_boundary_cell() -> None:
    """A relaxed WBM lattice on a Niggli boundary matches itself with zero RMSD."""
    lattice = [
        [3.5162634007, 1.7628e-6, -1.02026e-5],
        [-5.28446e-5, 8.4037724935, 4.8518887469],
        [-2.1294e-6, 3.8997e-6, 9.7038417713],
    ]
    structures = {
        "boundary": Structure(lattice, ["Si", "O"], [[0, 0, 0], [0.25, 0.26, 0.2]])
    }
    summary = symmetry.get_sym_info_from_structs(structures)
    result = symmetry.pred_vs_ref_struct_symmetry(
        summary, summary, structures, structures
    )
    # exactly 0 on macOS, 4.2e-17 on Linux CI (BLAS round-off, below f64 eps 2.2e-16);
    # a boundary mismatch gives None or a large RMSD
    rmsd = result.loc["boundary", Key.structure_rmsd_vs_dft]
    assert rmsd == pytest.approx(0, abs=1e-14)


def test_analyze_symmetry_failure_names_structure(cubic_struct: Structure) -> None:
    """A failed symmetry search raises with a note naming the material ID."""
    overlapping = Structure(Lattice.cubic(4), ["Na", "Cl"], np.zeros((2, 3)))
    with pytest.raises(ValueError, match="symmetry search failed") as exc_info:
        symmetry.get_sym_info_from_structs(
            {"ok": cubic_struct, "overlapping": overlapping}, pbar=False
        )
    assert exc_info.value.__notes__ == [f"{Key.mat_id}='overlapping'"]


@pytest.mark.parametrize("cell", ["supercell", "primitive"])
def test_analyze_symmetry_is_cell_choice_invariant(
    cubic_struct: Structure, cell: str
) -> None:
    """Supercells and primitive cells keep the conventional cell's space group."""
    other_struct = (
        cubic_struct * (2, 2, 2)
        if cell == "supercell"
        else SpacegroupAnalyzer(cubic_struct).get_primitive_standard_structure()
    )
    df_conventional = symmetry.get_sym_info_from_structs({"conventional": cubic_struct})
    df_other = symmetry.get_sym_info_from_structs({cell: other_struct})
    for df_sym, mat_id in ((df_conventional, "conventional"), (df_other, cell)):
        assert df_sym.index.name == Key.mat_id
        assert df_sym.index.tolist() == [mat_id]
        # space group name is a property of the lattice, not of the cell choice
        assert df_sym[Key.international_spg_name].iloc[0] == "I m -3 m"

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


def test_analyze_symmetry_with_ase_atoms(cubic_struct: Structure) -> None:
    """ASE Atoms give the right symmetry and match equivalent pymatgen Structures."""
    from ase import Atoms
    from ase.spacegroup import crystal

    atoms = crystal("Cu", [(0, 0, 0)], spacegroup=225, cellpar=[3.6] * 3 + [90] * 3)
    assert isinstance(atoms, Atoms)

    pbar_kwargs = {"disable": True}
    df_ase = symmetry.get_sym_info_from_structs({"ase": atoms}, pbar=pbar_kwargs)
    assert pbar_kwargs == {"disable": True}  # caller's dict is not mutated

    assert df_ase.index.tolist() == ["ase"]
    assert df_ase.index.name == Key.mat_id
    fcc_row = df_ase.loc["ase"]
    # conventional 4-atom fcc cell: 48 point ops x 4 centering translations
    assert fcc_row[Key.spg_num] == 225
    assert fcc_row[Key.international_spg_name] == "F m -3 m"
    assert fcc_row[Key.hall_symbol] == "-F 4 2 3"
    assert fcc_row[Key.n_sym_ops] == 192
    assert fcc_row[Key.wyckoff_symbols] == ["a"] * 4

    df_mixed = symmetry.get_sym_info_from_structs({"pmg": cubic_struct, "ase": atoms})
    assert df_mixed.index.tolist() == ["pmg", "ase"]
    assert df_mixed[Key.spg_num].tolist() == [229, 225]

    # should give same results for Structure and equivalent Atoms
    df_struct = symmetry.get_sym_info_from_structs({"struct": cubic_struct})
    df_atoms = symmetry.get_sym_info_from_structs(
        {"atoms": cubic_struct.to_ase_atoms()}
    )
    df_struct.index = df_atoms.index
    pd.testing.assert_frame_equal(df_struct, df_atoms)
