"""Functions to analyze symmetry of sets of structures."""

from collections.abc import Mapping
from typing import Any

import pandas as pd
from kemio import ItemError, Structure, rms_distances, symmetry_summary
from kemio.interop import to_structure
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey


def get_sym_info_from_structs(
    structures: Mapping[str, Structure | Mapping[str, Any]],
    *,
    symprec: float = 1e-2,
    angle_tolerance: float | None = None,
) -> pd.DataFrame:
    """Compile DataFrame of high-level symmetry information for a dictionary of
    structures, searching their symmetry in parallel (kemio.symmetry_summary). Rows
    whose symmetry search fails are NaN/None (calc_geo_opt_metrics skips them).

    Args:
        structures (dict[str, Structure | dict]): Map of material IDs to structures
            (or pymatgen-style structure dicts)
        symprec (float, optional): Symmetry precision in Å. Defaults to 1e-2.
        angle_tolerance (float, optional): Angle tolerance in radians (unlike spglib
            which uses degrees!). Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing symmetry information for each structure
    """
    mat_ids = list(structures)
    structs = []
    for mat_id, struct in structures.items():
        try:
            structs.append(to_structure(struct))
        except TypeError as exc:
            raise TypeError(f"{exc} for {mat_id!r}") from None
    try:
        summary = symmetry_summary(
            structs, symprec=symprec, angle_tolerance=angle_tolerance
        )
    except ItemError as exc:  # e.g. a structure without sites
        exc.add_note(f"{Key.mat_id}={mat_ids[exc.index]!r}")
        raise
    # nullable integers: a failed search leaves NA, not a float column (225.0 in CSVs)
    numbers, hall_numbers, n_ops = (
        pd.array(column, dtype="Int64")
        for column in (summary.number, summary.hall_number, summary.n_operations)
    )
    return pd.DataFrame(
        {
            Key.spg_num: numbers,
            Key.hall_num: hall_numbers,
            MbdKey.international_spg_name: summary.hm_short,
            Key.wyckoff_symbols: summary.wyckoffs,
            Key.n_sym_ops: n_ops,
            # every operation pairs one rotation with one translation
            Key.n_rot_syms: n_ops,
            Key.n_trans_syms: n_ops,
            Key.hall_symbol: summary.hall_symbol,
            Key.symprec: symprec,
            Key.angle_tolerance: angle_tolerance,
        },
        index=pd.Index(mat_ids, name=Key.mat_id),
    )


def pred_vs_ref_struct_symmetry(
    df_sym_pred: pd.DataFrame,
    df_sym_ref: pd.DataFrame,
    pred_structs: Mapping[str, Structure | Mapping[str, Any]],
    ref_structs: Mapping[str, Structure | Mapping[str, Any]],
) -> pd.DataFrame:
    """Get RMSD and compare symmetry between ML and DFT reference structures.

    Returns a copy of df_sym_pred with added columns for symmetry differences and
    RMSDs; the input DataFrame is left unchanged.

    Args:
        df_sym_pred (pd.DataFrame): symmetry information for ML model as returned by
            analyze_symmetry.
        df_sym_ref (pd.DataFrame): symmetry information for DFT reference as returned by
            analyze_symmetry.
        pred_structs (dict[str, Structure]): Map material IDs to ML-relaxed structures
        ref_structs (dict[str, Structure]): Map material IDs to reference structures

    Returns:
        pd.DataFrame: with added columns for symmetry differences
    """
    if df_sym_ref.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_ref.index.name=} must be {Key.mat_id!s}")
    if df_sym_pred.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_pred.index.name=} must be {Key.mat_id!s}")

    df_result = df_sym_pred.copy()

    df_result[MbdKey.spg_num_diff] = df_sym_pred[Key.spg_num] - df_sym_ref[Key.spg_num]
    df_result[MbdKey.n_sym_ops_diff] = (
        df_sym_pred[Key.n_sym_ops] - df_sym_ref[Key.n_sym_ops]
    )

    # scale=False and stol=1 are important for getting accurate distance of atomic
    # positions from DFT-relaxed positions. details in https://github.com/janosh/matbench-discovery/issues/230
    ref_ids, pred_ids = set(ref_structs), set(pred_structs)
    shared_ids = sorted(ref_ids & pred_ids)
    if len(shared_ids) == 0:
        raise ValueError(f"No shared IDs between:\n{pred_ids=}\n{ref_ids=}")

    pairs = [
        (to_structure(pred_structs[mat_id]), to_structure(ref_structs[mat_id]))
        for mat_id in shared_ids
    ]
    # one parallel call returning normalized RMS and max displacement (NaN for no
    # match). Individual values can differ from pymatgen StructureMatcher.
    nan_pair = (float("nan"), float("nan"))
    df_rms = pd.DataFrame(
        [match or nan_pair for match in rms_distances(pairs, stol=1.0)],
        index=shared_ids,
        columns=[MbdKey.structure_rmsd_vs_dft, Key.max_pair_dist],
        dtype=float,
    )
    df_result[df_rms.columns] = df_rms.reindex(df_result.index)

    return df_result
