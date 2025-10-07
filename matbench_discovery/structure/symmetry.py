"""Functions to analyze symmetry of sets of structures."""

from typing import Any

import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatviz.enums import Key
from pymatviz.typing import AnyStructure
from tqdm import tqdm

from matbench_discovery.enums import MbdKey


def get_sym_info_from_structs(
    structures: dict[str, AnyStructure],
    *,
    pbar: bool | dict[str, Any] = True,
    symprec: float = 1e-2,
    angle_tolerance: float | None = None,
) -> pd.DataFrame:
    """Compile DataFrame of high-level symmetry information for a dictionary of
    structures.

    Args:
        structures (dict[str, Structure | Atoms]): Map of material IDs to pymatgen
            Structures or ASE Atoms objects
        pbar (bool | dict[str, Any], optional): Whether to show progress bar. Defaults
            to True.
        symprec (float, optional): Symmetry precision of moyopy. Defaults to 1e-2.
        angle_tolerance (float, optional): Angle tolerance of moyopy (in radians unlike
            spglib which uses degrees!). Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing symmetry information for each structure
    """
    import moyopy
    from moyopy.interface import MoyoAdapter

    results: dict[str, dict[str, str | int | list[str] | float | None]] = {}
    iterator = structures.items()
    if pbar:
        pbar_kwargs = pbar if isinstance(pbar, dict) else {}
        pbar_kwargs.setdefault("desc", "Analyzing symmetry")
        iterator = tqdm(iterator, total=len(structures), **pbar_kwargs)

    for struct_key, struct in iterator:
        moyo_cell = MoyoAdapter.from_py_obj(struct)

        sym_data = moyopy.MoyoDataset(
            moyo_cell, symprec=symprec, angle_tolerance=angle_tolerance
        )

        sym_ops = sym_data.operations
        hall_symbol_entry = moyopy.HallSymbolEntry(hall_number=sym_data.hall_number)

        sym_info = {
            Key.spg_num: sym_data.number,
            Key.hall_num: sym_data.hall_number,
            MbdKey.international_spg_name: sym_data.site_symmetry_symbols,
            Key.wyckoff_symbols: sym_data.wyckoffs,
            Key.n_sym_ops: sym_ops.num_operations,
            Key.n_rot_syms: len(sym_ops.rotations),
            Key.n_trans_syms: len(sym_ops.translations),
            Key.hall_symbol: hall_symbol_entry.hm_short,
            Key.hall_num: hall_symbol_entry.hall_number,
        }
        results[struct_key] = sym_info | dict(
            symprec=symprec, angle_tolerance=angle_tolerance
        )

    df_sym = pd.DataFrame(results).T
    df_sym.index.name = Key.mat_id
    return df_sym


def pred_vs_ref_struct_symmetry(
    df_sym_pred: pd.DataFrame,
    df_sym_ref: pd.DataFrame,
    pred_structs: dict[str, Structure],
    ref_structs: dict[str, Structure],
    *,
    pbar: bool | dict[str, Any] = True,
) -> pd.DataFrame:
    """Get RMSD and compare symmetry between ML and DFT reference structures.

    Modifies the df_sym_pred DataFrame in place by adding columns for symmetry
    differences and RMSDs.

    Args:
        df_sym_pred (pd.DataFrame): symmetry information for ML model as returned by
            analyze_symmetry.
        df_sym_ref (pd.DataFrame): symmetry information for DFT reference as returned by
            analyze_symmetry.
        pred_structs (dict[str, Structure]): Map material IDs to ML-relaxed structures
        ref_structs (dict[str, Structure]): Map material IDs to reference structures
        pbar (bool | dict[str, Any], optional): Whether to show progress bar. Defaults
            to True.

    Returns:
        pd.DataFrame: with added columns for symmetry differences
    """
    if df_sym_ref.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_ref.index.name=} must be {Key.mat_id!s}")
    if df_sym_pred.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_pred.index.name=} must be {Key.mat_id!s}")

    df_result = df_sym_pred.copy()

    # Calculate differences
    df_result[MbdKey.spg_num_diff] = df_sym_pred[Key.spg_num] - df_sym_ref[Key.spg_num]
    df_result[MbdKey.n_sym_ops_diff] = (
        df_sym_pred[Key.n_sym_ops] - df_sym_ref[Key.n_sym_ops]
    )

    # scale=False and stol=1 are important for getting accurate distance of atomic
    # positions from DFT-relaxed positions. details in https://github.com/janosh/matbench-discovery/issues/230
    structure_matcher = StructureMatcher(stol=1.0, scale=False)
    ref_ids, pred_ids = set(ref_structs), set(pred_structs)
    shared_ids = ref_ids & pred_ids
    if len(shared_ids) == 0:
        raise ValueError(f"No shared IDs between:\n{pred_ids=}\n{ref_ids=}")

    # Initialize RMSD column
    df_result[MbdKey.structure_rmsd_vs_dft] = None

    if pbar:
        pbar_kwargs = pbar if isinstance(pbar, dict) else {}
        shared_ids = tqdm(
            shared_ids,
            **dict(leave=False, desc="Calculating RMSD") | pbar_kwargs,
        )

    for mat_id in shared_ids:
        rmsd, max_dist = structure_matcher.get_rms_dist(
            pred_structs[mat_id], ref_structs[mat_id]
        ) or (None, None)
        df_result.loc[mat_id, MbdKey.structure_rmsd_vs_dft] = rmsd
        df_result.loc[mat_id, Key.max_pair_dist] = max_dist

    return df_result
