"""Functions to calculate and save geometry optimization metrics."""

from typing import Any, Final

import numpy as np
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

# normalized RMSD at and beyond which a structure scores 0 in rmsd_auc
RMSD_AUC_CUTOFF: Final = 0.05


def rmsd_auc(
    rmsds: pd.Series, ref_ids: pd.Index, cutoff: float = RMSD_AUC_CUTOFF
) -> float:
    """Normalized area under the RMSD CDF on [0, cutoff] over all reference structures.

    The integral of the empirical CDF F over [0, c] equals mean(max(c - rmsd, 0)), so
    this is mean(max(1 - rmsd / cutoff, 0)), exact and independent of how a plotted CDF
    is sampled. Every reference structure counts: a missing prediction or a structure
    the matcher could not pair with its DFT counterpart scores 0, the same as one at
    or beyond the cutoff. Predictions for IDs outside ref_ids are not part of the
    test set and are ignored. 1 means every structure reproduces DFT exactly.
    """
    if len(ref_ids) == 0:
        raise ValueError("rmsd_auc needs at least one reference structure")
    rmsd = pd.to_numeric(rmsds.reindex(ref_ids), errors="raise").to_numpy(float)
    scores = np.clip(1 - rmsd / cutoff, 0, None)
    return float(np.nan_to_num(scores, nan=0.0).mean())


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

    # RMSD is symprec-invariant, so its score lives at task level, not per symprec
    score_key = str(MbdKey.rmsd_auc)
    task_scores = {score_key: round(float(df_geo_opt[score_key].iloc[0]), 4)}
    update_yaml_file(
        model.yaml_path,
        "metrics.geo_opt",
        commented_map_with_units(task_scores, {score_key: "fraction"}),
    )

    metrics_for_symprec = {
        str(output_key): round(float(df_geo_opt[input_key].iloc[0]), 4)
        for output_key, input_key in (
            (Key.n_sym_ops_mae, Key.n_sym_ops_mae),
            (Key.symmetry_decrease, Key.symmetry_decrease),
            (Key.symmetry_match, Key.symmetry_match),
            (Key.symmetry_increase, Key.symmetry_increase),
            (MbdKey.symmetry_lateral, MbdKey.symmetry_lateral),
        )
    }
    metrics_for_symprec[str(Key.n_structures)] = int(
        df_geo_opt[Key.n_structures].iloc[0]
    )
    metric_units: dict[str, str] = {
        Key.n_sym_ops_mae: "unitless",
        Key.symmetry_decrease: "fraction",
        Key.symmetry_match: "fraction",
        Key.symmetry_increase: "fraction",
        MbdKey.symmetry_lateral: "fraction",
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


def calc_geo_opt_metrics(
    df_model_analysis: pd.DataFrame, ref_ids: pd.Index
) -> dict[str, float]:
    """Calculate geometry optimization metrics for a single model.

    Args:
        df_model_analysis (pd.DataFrame): DataFrame with geometry optimization metrics
            for one model, indexed by material ID. Required columns are:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT
        ref_ids (pd.Index): Material IDs of the DFT reference (the test set). Rows
            outside it are dropped; reference IDs the model lacks count as failures.

    Returns:
        dict[str, float]: Geometry optimization metrics with keys:
            - rmsd_auc: Normalized RMSD CDF area up to RMSD_AUC_CUTOFF (see rmsd_auc)
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease: Fraction of structures with decreased symmetry
            - symmetry_match: Fraction of structures with matching symmetry
            - symmetry_increase: Fraction of structures with increased symmetry
            - symmetry_lateral: Fraction whose space group changed with the same
              number of symmetry operations (the 4 fractions sum to 1)
            - n_structures: Number of structures evaluated

    Notes:
        - total number of structures is counted based on valid symmetry data
        - symmetry metrics are calculated only on structures with valid symmetry data
    """
    df_model_analysis = df_model_analysis.reindex(ref_ids)
    spg_diff = df_model_analysis[MbdKey.spg_num_diff]
    n_sym_ops_diff = df_model_analysis[MbdKey.n_sym_ops_diff]

    # symmetry metrics use only structures with valid symmetry results: detection can
    # fail in the symmetry finder itself, leaving model vs algo blame unassignable
    valid_sym_mask = spg_diff.notna()
    n_valid_sym = valid_sym_mask.sum()

    sym_ops_mae = n_sym_ops_diff[valid_sym_mask].abs().mean()

    changed_mask = (spg_diff != 0) & valid_sym_mask
    sym_decreased = (n_sym_ops_diff < 0) & changed_mask
    sym_increased = (n_sym_ops_diff > 0) & changed_mask
    sym_lateral = (n_sym_ops_diff == 0) & changed_mask
    sym_matched = ~changed_mask & valid_sym_mask

    return {
        str(MbdKey.rmsd_auc): rmsd_auc(
            df_model_analysis[MbdKey.structure_rmsd_vs_dft], ref_ids
        ),
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
        str(MbdKey.symmetry_lateral): float(sym_lateral.sum() / n_valid_sym)
        if n_valid_sym > 0
        else float("nan"),
        str(Key.n_structures): n_valid_sym,
    }
