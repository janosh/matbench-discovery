"""Functions to calculate and save geometry optimization metrics."""

from typing import Any

import pandas as pd
from pymatviz.enums import Key

from matbench_discovery import repo_relative_path
from matbench_discovery.data import (
    canonical_scientific_notation,
    commented_map_with_units,
    merge_file_ref,
    update_yaml_file,
)
from matbench_discovery.enums import MbdKey, Model


def write_metrics_to_yaml(
    df_geo_opt: pd.DataFrame,
    model: Model,
    symprec: float,
    analysis_file_path: str,
) -> dict[str, object]:
    """Write geometric optimization metrics to model's YAML file.

    The analysis_file ref keeps its Figshare url/size/md5 while the file name is
    unchanged (e.g. metrics re-derived from a cached CSV); a renamed analysis file
    (new date or moyo version) drops them until the next upload.

    Args:
        df_geo_opt (pd.DataFrame): Geometric optimization metrics
        model (Model): Model to write metrics for
        symprec (float): Symmetry precision used for analysis
        analysis_file_path (str): Path to analysis file

    Returns:
        dict[str, object]: Geometric optimization metrics for this model and
            symmetry precision.
    """
    analysis_file_path = repo_relative_path(analysis_file_path)

    metrics_for_symprec = {
        str(Key.rmsd): round(
            float(df_geo_opt[MbdKey.structure_rmsd_vs_dft].iloc[0]), 4
        ),
        str(Key.n_sym_ops_mae): round(float(df_geo_opt[Key.n_sym_ops_mae].iloc[0]), 4),
        str(Key.symmetry_decrease): round(
            float(df_geo_opt[Key.symmetry_decrease].iloc[0]), 4
        ),
        str(Key.symmetry_match): round(
            float(df_geo_opt[Key.symmetry_match].iloc[0]), 4
        ),
        str(Key.symmetry_increase): round(
            float(df_geo_opt[Key.symmetry_increase].iloc[0]), 4
        ),
        str(Key.n_structures): int(df_geo_opt[Key.n_structures].iloc[0]),
    }
    metric_units: dict[str, str] = {
        Key.rmsd: "unitless",
        Key.n_sym_ops_mae: "unitless",
        Key.symmetry_decrease: "fraction",
        Key.symmetry_match: "fraction",
        Key.symmetry_increase: "fraction",
        Key.n_structures: "count",
    }
    metrics_for_symprec = commented_map_with_units(metrics_for_symprec, metric_units)
    symprec_key = f"symprec={canonical_scientific_notation(symprec)}"

    def merge_block(previous: dict[str, Any]) -> dict[str, Any]:
        """Merge new metrics over the prior block under update_yaml_file's lock."""
        metrics_for_symprec["analysis_file"] = merge_file_ref(
            previous.get("analysis_file"), analysis_file_path
        )
        for key, val in previous.items():  # keep prior keys not emitted here
            metrics_for_symprec.setdefault(key, val)
        return metrics_for_symprec

    update_yaml_file(model.yaml_path, f"metrics.geo_opt.{symprec_key}", merge_block)
    return metrics_for_symprec


def calc_geo_opt_metrics(df_model_analysis: pd.DataFrame) -> dict[str, float]:
    """Calculate geometry optimization metrics for a single model.

    Args:
        df_model_analysis (pd.DataFrame): DataFrame with geometry optimization metrics
            for one model. Required columns are:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT

    Returns:
        dict[str, float]: Geometry optimization metrics with keys:
            - structure_rmsd_vs_dft: Mean RMSD between predicted and DFT structures
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease: Fraction of structures with decreased symmetry
            - symmetry_match: Fraction of structures with matching symmetry
            - symmetry_increase: Fraction of structures with increased symmetry
            - n_structures: Number of structures evaluated

    Notes:
        - total number of structures is counted based on valid symmetry data
        - NaN RMSD values are filled with 1.0 (the stol value set in StructureMatcher)
        - symmetry metrics are calculated only on structures with valid symmetry data
    """
    spg_diff = df_model_analysis[MbdKey.spg_num_diff]
    n_sym_ops_diff = df_model_analysis[MbdKey.n_sym_ops_diff]
    rmsd_vals = df_model_analysis[MbdKey.structure_rmsd_vs_dft]

    # symmetry metrics use only structures with valid symmetry results: detection can
    # fail in the symmetry finder itself, leaving model vs algo blame unassignable
    valid_sym_mask = spg_diff.notna()
    n_valid_sym = valid_sym_mask.sum()

    # Fill NaN values with 1.0 (the stol value we set in StructureMatcher)
    mean_rmsd = pd.to_numeric(rmsd_vals, errors="coerce").fillna(1.0).mean()

    sym_ops_mae = n_sym_ops_diff[valid_sym_mask].abs().mean()

    changed_mask = (spg_diff != 0) & valid_sym_mask
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
