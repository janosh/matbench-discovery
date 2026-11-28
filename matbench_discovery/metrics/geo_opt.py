"""Functions to calculate and save geometry optimization metrics."""

from typing import Any, Final

import numpy as np
import pandas as pd

from matbench_discovery import repo_relative_path
from matbench_discovery.data import (
    canonical_scientific_notation,
    commented_map_with_units,
    merge_file_ref,
    update_yaml_file,
)
from matbench_discovery.enums import Key, Model

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
    test set and are ignored. 1 means every structure's atomic positions match DFT
    exactly. RMSD is blind to lattice strain the matcher accepts (ltol 0.2, angle_tol
    5 degrees): such a cell scores as if it matched DFT's.

    Raises:
        ValueError: On an empty or duplicated test set, duplicate prediction IDs, no
            overlap between prediction and test-set IDs (e.g. a CSV read without
            index_col), non-float RMSDs, or negative RMSDs.
    """
    if len(ref_ids) == 0:
        raise ValueError("rmsd_auc needs at least one reference structure")
    for label, ids in (("ref_ids", ref_ids), ("rmsds.index", rmsds.index)):
        if dupes := ids[ids.duplicated()].unique().tolist():
            raise ValueError(f"duplicate {label} {dupes[:5]}")
    if not rmsds.index.isin(ref_ids).any():
        raise ValueError(
            f"no RMSD index entries are reference IDs: {rmsds.index[:3].tolist()} vs "
            f"{ref_ids[:3].tolist()}"
        )
    if not pd.api.types.is_float_dtype(rmsds.dtype):
        raise ValueError(f"RMSDs must be floats (NaN = unmatched), got {rmsds.dtype}")
    rmsd = rmsds.reindex(ref_ids).to_numpy(float)
    if (rmsd < 0).any():  # NaN compares False
        raise ValueError(f"negative RMSDs: {rmsd[rmsd < 0][:5].tolist()}")
    scores = np.clip(1 - rmsd / cutoff, 0, None)
    return float(np.nan_to_num(scores, nan=0.0).mean())


def write_metrics_to_yaml(
    df_geo_opt: pd.DataFrame,
    model: Model,
    symprec: float,
    analysis_file_path: str,
    *,
    regenerated: bool = False,
) -> dict[str, object]:
    """Write geometric optimization metrics to model's YAML file.

    A cached file retains its registered URL and checksum. A renamed or regenerated
    CSV gets a fresh file reference until its bytes are published.

    Args:
        df_geo_opt (pd.DataFrame): Geometric optimization metrics
        model (Model): Model to write metrics for
        symprec (float): Symmetry precision used for analysis
        analysis_file_path (str): Path to analysis file
        regenerated (bool): Whether this run wrote new analysis bytes.

    Returns:
        dict[str, object]: Geometric optimization metrics for this model and
            symmetry precision.
    """
    analysis_file_path = repo_relative_path(analysis_file_path)

    row = df_geo_opt.iloc[0]
    # RMSD is symprec-invariant, so its score lives at task level, not per symprec
    score_key = str(Key.rmsd_auc)
    task_scores = {score_key: round(float(row[score_key]), 4)}
    update_yaml_file(
        model.yaml_path,
        "metrics.geo_opt",
        commented_map_with_units(task_scores, {score_key: "fraction"}),
    )

    metric_units: dict[str, str] = {
        Key.n_sym_ops_mae: "unitless",
        Key.symmetry_decrease: "fraction",
        Key.symmetry_match: "fraction",
        Key.symmetry_increase: "fraction",
        Key.symmetry_lateral: "fraction",
        Key.n_structures: "count",
    }
    metrics_for_symprec = commented_map_with_units(
        {
            str(key): int(row[key]) if unit == "count" else round(float(row[key]), 4)
            for key, unit in metric_units.items()
        },
        metric_units,
    )
    symprec_key = f"symprec={canonical_scientific_notation(symprec)}"

    def merge_block(previous: dict[str, Any]) -> dict[str, Any]:
        """Merge new metrics over the prior block under update_yaml_file's lock."""
        metrics_for_symprec["analysis_file"] = merge_file_ref(
            previous.get("analysis_file"), analysis_file_path, replace=regenerated
        )
        for key, val in previous.items():  # keep prior keys not emitted here
            metrics_for_symprec.setdefault(key, val)
        return metrics_for_symprec

    update_yaml_file(model.yaml_path, f"metrics.geo_opt.{symprec_key}", merge_block)
    model.__dict__.pop("metadata", None)
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
            - n_structures: Number of structures with valid symmetry data, the
              only ones the symmetry metrics are calculated on
    """
    df_model_analysis = df_model_analysis.reindex(ref_ids)
    spg_diff = df_model_analysis[Key.spg_num_diff]
    n_sym_ops_diff = df_model_analysis[Key.n_sym_ops_diff]

    # symmetry metrics use only structures with valid symmetry results: detection can
    # fail in the symmetry finder itself, leaving model vs algo blame unassignable
    valid_sym_mask = spg_diff.notna()
    n_valid_sym = valid_sym_mask.sum()
    changed_mask = (spg_diff != 0) & valid_sym_mask
    sym_buckets = {
        Key.symmetry_decrease: (n_sym_ops_diff < 0) & changed_mask,
        Key.symmetry_match: ~changed_mask & valid_sym_mask,
        Key.symmetry_increase: (n_sym_ops_diff > 0) & changed_mask,
        Key.symmetry_lateral: (n_sym_ops_diff == 0) & changed_mask,
    }

    return {
        str(Key.rmsd_auc): rmsd_auc(
            df_model_analysis[Key.structure_rmsd_vs_dft], ref_ids
        ),
        str(Key.n_sym_ops_mae): float(n_sym_ops_diff[valid_sym_mask].abs().mean()),
        **{
            str(key): float(mask.sum() / n_valid_sym) if n_valid_sym else float("nan")
            for key, mask in sym_buckets.items()
        },
        str(Key.n_structures): n_valid_sym,
    }
