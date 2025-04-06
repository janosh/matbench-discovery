# ruff: noqa: FBT001, FBT002
import argparse
import importlib
import importlib.metadata
import multiprocessing as mp
import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from typing import Final

import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.data import DataFiles


def get_sym_info_from_structs(
    structures: dict[str, Structure],
    *,
    pbar: bool | dict[str, str | float | bool] = True,
    symprec: float = 1e-2,
    angle_tolerance: float | None = None,
) -> pd.DataFrame:
    """Compile DataFrame of high-level symmetry information for a dictionary of
    structures.

    Args:
        structures (dict[str, Structure | Atoms]): Map of material IDs to pymatgen
            Structures or ASE Atoms objects
        pbar (bool | dict[str, str | float | bool], optional): Whether to show progress
            bar. Defaults to True.
        symprec (float, optional): Symmetry precision of moyopy. Defaults to 1e-2.
        angle_tolerance (float, optional): Angle tolerance of moyopy (in radians unlike
            spglib which uses degrees!). Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing symmetry information for each structure
    """
    import moyopy
    from moyopy.interface import MoyoAdapter

    results: dict[str, dict[str, str | int | list[str]]] = {}
    iterator = structures.items()
    if pbar:
        pbar_kwargs = pbar if isinstance(pbar, dict) else {}
        pbar_kwargs.setdefault("desc", "Analyzing symmetry")
        iterator = tqdm(iterator, total=len(structures), **pbar_kwargs)

    for struct_key, struct in iterator:
        moyo_cell = MoyoAdapter.from_py_obj(struct)

        try:
            sym_data = moyopy.MoyoDataset(
                moyo_cell, symprec=symprec, angle_tolerance=angle_tolerance
            )
            sym_ops = sym_data.operations
            hall_symbol_entry = moyopy.HallSymbolEntry(hall_number=sym_data.hall_number)

            sym_info = {
                "spg_num": sym_data.number,
                "international_spg_name": sym_data.site_symmetry_symbols,
                "wyckoff_symbols": sym_data.wyckoffs,
                "n_sym_ops": sym_ops.num_operations,
                "n_rot_syms": len(sym_ops.rotations),
                "n_trans_syms": len(sym_ops.translations),
                "hall_symbol": hall_symbol_entry.hm_short,
                "hall_num": hall_symbol_entry.hall_number,
            }

        except Exception as e:
            print(f"Error analyzing symmetry for {struct_key}: {e}")
            sym_info = None

        if sym_info is None:
            sym_info = dict(symprec=symprec, angle_tolerance=angle_tolerance)
        results[struct_key] = sym_info

    df_sym = pd.DataFrame(results).T
    df_sym.index.name = Key.mat_id
    return df_sym


def process_structure(
    id_pred_ref_tuple: tuple[str, Structure, Structure],
) -> tuple[str, float, float]:
    mat_id, pred_struct, ref_struct = id_pred_ref_tuple
    # scale=False and stol=1 are important for getting accurate distance of atomic
    # positions from DFT-relaxed positions. details in https://github.com/janosh/matbench-discovery/issues/230
    structure_matcher = StructureMatcher(stol=1.0, scale=False)
    rmsd, max_dist = structure_matcher.get_rms_dist(pred_struct, ref_struct) or (
        None,
        None,
    )
    return mat_id, rmsd, max_dist  # type: ignore[return-value]


def pred_vs_ref_struct_symmetry(
    df_sym_pred: pd.DataFrame,
    df_sym_ref: pd.DataFrame,
    pred_structs: dict[str, Structure],
    ref_structs: dict[str, Structure],
    *,
    pbar: bool | dict[str, str | float | bool] = True,
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
        pbar (bool | dict[str, str | float | bool], optional): Whether to show progress
            bar. Defaults to True.

    Returns:
        pd.DataFrame: with added columns for symmetry differences
    """
    if df_sym_ref.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_ref.index.name=} must be {Key.mat_id!s}")
    if df_sym_pred.index.name != Key.mat_id:
        raise ValueError(f"{df_sym_pred.index.name=} must be {Key.mat_id!s}")

    df_result = df_sym_pred.copy()

    # Calculate differences
    df_result["spg_num_diff"] = df_sym_pred["spg_num"] - df_sym_ref["spg_num"]
    df_result["n_sym_ops_diff"] = df_sym_pred["n_sym_ops"] - df_sym_ref["n_sym_ops"]

    ref_ids, pred_ids = set(ref_structs), set(pred_structs)
    shared_ids = ref_ids & pred_ids
    if len(shared_ids) == 0:
        raise ValueError(f"No shared IDs between:\n{pred_ids=}\n{ref_ids=}")

    # Initialize RMSD column
    df_result["structure_rmsd_vs_dft"] = None

    id_pred_ref_tuples = [
        (mat_id, pred_structs[mat_id], ref_structs[mat_id]) for mat_id in shared_ids
    ]

    if pbar:
        pbar_kwargs = pbar if isinstance(pbar, dict) else {}
        id_pred_ref_tuples = tqdm(
            id_pred_ref_tuples, **dict(desc="Calculating RMSD") | pbar_kwargs
        )

    # Process structures in parallel with 10 workers
    with ProcessPoolExecutor(max_workers=10) as executor:
        results = executor.map(process_structure, id_pred_ref_tuples)

    # Update the DataFrame with results
    for mat_id, rmsd, max_dist in results:
        df_result.loc[mat_id, "structure_rmsd_vs_dft"] = rmsd
        df_result.loc[mat_id, "max_pair_dist"] = max_dist

    return df_result


def calc_geo_opt_metrics(df_model_analysis: pd.DataFrame) -> dict[str, float]:
    """Calculate geometry optimization metrics for a single model.

    Args:
        df_model_analysis (pd.DataFrame): DataFrame with geometry optimization metrics
            for one model. Required columns are:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT
        model_name (str): Name of the model being analyzed.

    Returns:
        dict[str, float]: Geometry optimization metrics with keys:
            - structure_rmsd_vs_dft: Mean RMSD between predicted and DFT structures
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease: Fraction of structures with decreased symmetry
            - symmetry_match: Fraction of structures with matching symmetry
            - symmetry_increase: Fraction of structures with increased symmetry
            - n_structs: Number of structures evaluated
    """
    # Get relevant columns
    spg_diff = df_model_analysis["spg_num_diff"]
    n_sym_ops_diff = df_model_analysis["n_sym_ops_diff"]
    rmsd_vals = df_model_analysis["structure_rmsd_vs_dft"]

    # Count total number of structures (excluding NaN values)
    n_structs = len(spg_diff.dropna())

    # Fill NaN values with 1.0 (the stol value we set in StructureMatcher)
    mean_rmsd = rmsd_vals.fillna(1.0).mean()
    sym_ops_mae = n_sym_ops_diff.abs().mean()

    # Count cases where spacegroup changed
    changed_mask = spg_diff != 0
    # Among changed cases, count whether symmetry increased or decreased
    sym_decreased = (n_sym_ops_diff < 0) & changed_mask
    sym_increased = (n_sym_ops_diff > 0) & changed_mask
    sym_matched = ~changed_mask

    return {
        "structure_rmsd_vs_dft": float(mean_rmsd),
        "n_sym_ops_mae": float(sym_ops_mae),
        "symmetry_decrease": float(sym_decreased.sum() / n_structs),
        "symmetry_match": float(sym_matched.sum() / n_structs),
        "symmetry_increase": float(sym_increased.sum() / n_structs),
        "n_structures": n_structs,
    }


def analyze_model_symprec(
    predicted_structures_path: str,
    symprec: float,
    moyo_version: str,
    df_dft_analysis: pd.DataFrame,
    dft_structs: dict[str, Structure],
    debug_mode: int = 0,
    pbar_pos: int = 0,  # tqdm progress bar position
    overwrite: bool = False,
    struct_col: str = "orb_structure",
) -> None:
    """Analyze a single model for a single symprec value."""

    # Load model structures
    if predicted_structures_path.endswith(
        (".json", ".jsonl", ".jsonl.gz", ".json.gz", ".json.xz")
    ):
        df_ml_structs = pd.read_json(predicted_structures_path)
        print("Columns:", df_ml_structs.columns)
    else:
        raise ValueError(
            "Relaxed structure analysis currently only supports pymatgen JSON, "
            f"got {predicted_structures_path}"
        )

    # try normalize material ID column or raise
    if Key.mat_id in df_ml_structs:
        df_ml_structs = df_ml_structs.set_index(Key.mat_id)
    elif df_ml_structs.index[0].startswith("wbm-"):
        df_ml_structs.index.name = Key.mat_id
    else:
        raise ValueError(f"Could not infer ID column from {df_ml_structs.columns}")

    if debug_mode:
        df_ml_structs = df_ml_structs.head(debug_mode)

    if struct_col not in df_ml_structs:
        struct_cols = [col for col in df_ml_structs if Key.structure in col]
        print(
            f"⚠️ {struct_col=} not found in structures loaded "
            f"from {predicted_structures_path}. Did you mean one of {struct_cols}?"
        )
        return

    # Convert structures
    model_structs = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_ml_structs[struct_col].items()
    }

    symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")
    geo_opt_filename = predicted_structures_path.removesuffix(".json.gz").removesuffix(
        ".jsonl.gz"
    )
    geo_opt_csv_path = f"{geo_opt_filename}-{symprec_str}-{moyo_version}.csv.gz"

    if os.path.isfile(geo_opt_csv_path):
        if overwrite:
            os.remove(geo_opt_csv_path)
        else:
            raise FileExistsError(
                f"{predicted_structures_path} already analyzed at {geo_opt_csv_path}"
            )

    # Analyze symmetry for current symprec
    pbar_desc = f"Process {pbar_pos}: Analyzing for {symprec=}"
    df_model_analysis = get_sym_info_from_structs(
        model_structs,
        pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
        symprec=symprec,
    )

    # Compare with DFT reference
    pbar_desc = f"Process {pbar_pos}:Comparing DFT vs predicted for {symprec=}"
    # break here
    df_ml_geo_analysis = pred_vs_ref_struct_symmetry(
        df_model_analysis,
        df_dft_analysis,
        model_structs,
        dft_structs,
        pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
    )

    # Save model results
    df_ml_geo_analysis.to_csv(geo_opt_csv_path)
    print(f"Complete and saved results to {geo_opt_csv_path}")

    # Calculate metrics and write to YAML
    metrics = calc_geo_opt_metrics(df_ml_geo_analysis)
    print(f"\nCalculated metrics: {metrics}")

    df_metrics = pd.DataFrame(metrics, index=[0])
    print(df_metrics.to_csv())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--predicted_structures_path",
        type=str,
        help="Path to predicted structures.",
    )
    parser.add_argument(
        "--symprec",
        nargs="+",
        type=float,
        default=[1e-2, 1e-5],
        help="Symmetry precision values to analyze.",
    )
    parser.add_argument(
        "--debug",
        type=int,
        default=0,
        help="If > 0, only analyze this many structures.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, mp.cpu_count() - 1),
        help="Number of processes to use. Defaults to number of model-symprec combos.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing analysis.",
    )
    parser.add_argument(
        "--struct_col",
        type=str,
        default="orb_structure",
        help="Column name in predicted structures file to use for structures.",
    )
    args, _unknown = parser.parse_known_args()

    predicted_structures_path = args.predicted_structures_path
    # set to > 0 to activate debug mode, only that many structures will be analyzed
    debug_mode: Final[int] = args.debug
    # List of symprec values to analyze
    symprec_values: Final[Sequence[float]] = args.symprec

    # Get list of models to analyze
    moyo_version = f"moyo={importlib.metadata.version('moyopy')}"

    print("Loading WBM PBE structures...")
    wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
    df_wbm_structs: pd.DataFrame = pd.read_json(wbm_cse_path).set_index(Key.mat_id)

    if debug_mode:
        df_wbm_structs = df_wbm_structs.head(debug_mode)

    dft_structs: dict[str, Structure] = {
        mat_id: Structure.from_dict(cse[Key.structure])
        for mat_id, cse in df_wbm_structs["computed_structure_entry"].items()
    }

    # Process DFT structures for each symprec value
    dft_analysis_dict: dict[float, pd.DataFrame] = {}
    for symprec in symprec_values:
        symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")
        dft_csv_path = (
            f"{ROOT}/data/wbm/dft-geo-opt-{symprec_str}-{moyo_version}.csv.gz"
        )

        if os.path.isfile(dft_csv_path):
            print(f"Loading DFT symmetries from {dft_csv_path}")
            dft_analysis_dict[symprec] = pd.read_csv(dft_csv_path).set_index(Key.mat_id)
            if len(dft_analysis_dict[symprec]) != len(dft_structs):
                raise ValueError(
                    f"Length mismatch in dft vs ref structures. "
                    f"{len(dft_analysis_dict[symprec])} != {len(dft_structs)}"
                )
        else:
            dft_analysis_dict[symprec] = get_sym_info_from_structs(
                dft_structs,
                pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
                symprec=symprec,
            )
            dft_analysis_dict[symprec].to_csv(dft_csv_path)

    # Create list of all model-symprec combinations
    n_workers = min(len(symprec_values), args.workers)

    # Process model-symprec combinations in parallel
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(
                analyze_model_symprec,
                predicted_structures_path,
                symprec,
                moyo_version,
                dft_analysis_dict[symprec],
                dft_structs,
                debug_mode,
                pbar_pos=idx,  # assign unique position to each task's progress bar
                overwrite=args.overwrite,
                struct_col=args.struct_col,
            )
            for idx, symprec in enumerate(symprec_values)
        ]
        # Wait for all tasks to complete
        for future in futures:
            future.result()  # This will raise any exceptions that occurred

    print(f"\nAll {len(symprec_values)} symprec values processed!")
