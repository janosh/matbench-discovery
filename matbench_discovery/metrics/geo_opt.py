"""Functions to calculate and save geometry optimization metrics."""

import pandas as pd
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery import ROOT
from matbench_discovery.data import update_yaml_file
from matbench_discovery.enums import MbdKey, Model


def write_metrics_to_yaml(
    df_geo_opt: pd.DataFrame,
    model: Model,
    symprec: float,
    analysis_file_path: str,
) -> dict[str, str | float]:
    """Write geometric optimization metrics to model's YAML file.

    Args:
        df_geo_opt (pd.DataFrame): Geometric optimization metrics
        model (Model): Model to write metrics for
        symprec (float): Symmetry precision used for analysis
        analysis_file_path (str): Path to analysis file

    Returns:
        dict[str, str | float]: Geometric optimization metrics for this model and
            symmetry precision.
    """
    # Convert absolute path to relative path if needed
    analysis_file_path = analysis_file_path.removeprefix(f"{ROOT}/")

    # Get metrics for this model
    metrics_for_symprec = {
        str(Key.rmsd): float(round(df_geo_opt[MbdKey.structure_rmsd_vs_dft], 4)),
        str(Key.n_sym_ops_mae): float(round(df_geo_opt[Key.n_sym_ops_mae], 4)),
        str(Key.symmetry_decrease): float(round(df_geo_opt[Key.symmetry_decrease], 4)),
        str(Key.symmetry_match): float(round(df_geo_opt[Key.symmetry_match], 4)),
        str(Key.symmetry_increase): float(round(df_geo_opt[Key.symmetry_increase], 4)),
        str(Key.n_structures): int(df_geo_opt[Key.n_structures]),
        "analysis_file": analysis_file_path,
        "analysis_file_url": None,  # to be filled after uploading to figshare
    }
    metrics_for_symprec = CommentedMap(metrics_for_symprec)
    symprec_key = f"{symprec=:.0e}".replace("e-0", "e-")

    # Define units for metrics
    metric_units = {
        Key.rmsd: "unitless",
        Key.n_sym_ops_mae: "unitless",
        Key.symmetry_decrease: "fraction",
        Key.symmetry_match: "fraction",
        Key.symmetry_increase: "fraction",
        Key.n_structures: "count",
    }
    # Add units as YAML end-of-line comments
    for key in metrics_for_symprec:
        if unit := metric_units.get(key):
            metrics_for_symprec.yaml_add_eol_comment(unit, key, column=1)

    update_yaml_file(
        model.yaml_path, f"metrics.geo_opt.{symprec_key}", metrics_for_symprec
    )
    return metrics_for_symprec


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

    Notes:
        - total number of structures (n_structs) is counted based on valid symmetry data
        - NaN RMSD values are filled with 1.0 (the stol value set in StructureMatcher)
        - symmetry metrics are calculated only on structures with valid symmetry data
    """
    # Get relevant columns
    spg_diff = df_model_analysis[MbdKey.spg_num_diff]
    n_sym_ops_diff = df_model_analysis[MbdKey.n_sym_ops_diff]
    rmsd_vals = df_model_analysis[MbdKey.structure_rmsd_vs_dft]

    # For symmetry metrics, we only use structures with valid symmetry results
    # in rare cases, symmetry detection may fail because of the symmetry finder
    # algorithm rather than something being wrong with the model-relaxed structure so
    # not clear how to assign blame for missing results between model and symmetry algo
    valid_sym_mask = spg_diff.notna()
    n_valid_sym = valid_sym_mask.sum()

    # Fill NaN values with 1.0 (the stol value we set in StructureMatcher)
    mean_rmsd = rmsd_vals.infer_objects(copy=False).fillna(1.0).mean()

    # Calculate symmetry metrics only on valid symmetry data
    sym_ops_mae = n_sym_ops_diff[valid_sym_mask].abs().mean()

    # Count cases where spacegroup changed
    changed_mask = (spg_diff != 0) & valid_sym_mask
    # Among changed cases, count whether symmetry increased or decreased
    sym_decreased = (n_sym_ops_diff < 0) & changed_mask
    sym_increased = (n_sym_ops_diff > 0) & changed_mask
    sym_matched = ~changed_mask & valid_sym_mask

    return {
        str(MbdKey.structure_rmsd_vs_dft): float(mean_rmsd),
        str(Key.n_sym_ops_mae): float(sym_ops_mae),
        str(Key.symmetry_decrease): float(sym_decreased.sum() / n_valid_sym)
        if n_valid_sym > 0
        else float("nan"),
        str(Key.symmetry_match): float(sym_matched.sum() / n_valid_sym)
        if n_valid_sym > 0
        else float("nan"),
        str(Key.symmetry_increase): float(sym_increased.sum() / n_valid_sym)
        if n_valid_sym > 0
        else float("nan"),
        str(Key.n_structures): n_valid_sym,
    }
