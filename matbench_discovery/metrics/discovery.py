"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.data import MAX_E_FORM_ERROR_THRESHOLD
from matbench_discovery.enums import MbdKey, Model, TestSubset


def classify_stable(
    each_true: Sequence[float | None] | pd.Series | np.ndarray,
    each_pred: Sequence[float | None] | pd.Series | np.ndarray,
    *,
    stability_threshold: float = STABILITY_THRESHOLD,
    fillna: bool = True,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Classify model stability predictions as true/false positive/negatives (usually
    w.r.t DFT-ground truth labels). All energies are assumed to be in eV/atom
    (but shouldn't really matter as long as they're consistent).

    Args:
        each_true (Sequence[float] | pd.Series): Ground truth energy above convex hull
            values.
        each_pred (Sequence[float] | pd.Series): Model-predicted energy above convex
            hull values.
        stability_threshold (float, optional): Maximum energy above convex hull
            for a material to still be considered stable. Usually 0, 0.05 or 0.1.
            Defaults to STABILITY_THRESHOLD, meaning a material has to be directly on
            the hull to be called stable. Negative values mean a material has to pull
            the known hull down by that amount to count as stable. Few materials lie
            below the known hull, so only negative values very close to 0 make sense.
        fillna (bool): Whether to fill NaNs as the model predicting unstable. Defaults
            to True.

    Returns:
        tuple[TP, FN, FP, TN]: Indices as pd.Series for true positives,
            false negatives, false positives and true negatives (in this order).

    Raises:
        ValueError: If sum of positive + negative preds doesn't add up to the total.
    """
    if len(each_true) != len(each_pred):
        raise ValueError(f"{len(each_true)=} != {len(each_pred)=}")

    each_true_arr = pd.to_numeric(pd.Series(each_true), errors="coerce")
    each_pred_arr = pd.to_numeric(pd.Series(each_pred), errors="coerce")

    if stability_threshold is None or np.isnan(stability_threshold):
        raise ValueError("stability_threshold must be a real number")
    actual_pos = each_true_arr <= stability_threshold
    actual_neg = each_true_arr > stability_threshold

    model_pos = each_pred_arr <= stability_threshold
    model_neg = each_pred_arr > stability_threshold

    if fillna:
        nan_mask = each_pred_arr.isna()
        # for in both the model's stable and unstable preds, fill NaNs as unstable
        model_pos[nan_mask] = False
        model_neg[nan_mask] = True

        n_pos, n_neg, total = model_pos.sum(), model_neg.sum(), len(each_pred)
        if n_pos + n_neg != total:
            raise ValueError(
                f"after filling NaNs, the sum of positive ({n_pos}) and negative "
                f"({n_neg}) predictions should add up to {total=}"
            )

    true_pos = actual_pos & model_pos
    false_neg = actual_pos & model_neg
    false_pos = actual_neg & model_pos
    true_neg = actual_neg & model_neg

    return true_pos, false_neg, false_pos, true_neg


def stable_metrics(
    each_true: Sequence[float | None] | pd.Series | np.ndarray,
    each_pred: Sequence[float | None] | pd.Series | np.ndarray,
    *,
    stability_threshold: float = STABILITY_THRESHOLD,
    fillna: bool = True,
) -> dict[str, float]:
    """Get a dictionary of stability prediction metrics. Mostly binary classification
    metrics, but also MAE, RMSE and R2.

    Args:
        each_true (Sequence[float | None] | pd.Series): true energy above convex hull
        each_pred (Sequence[float | None] | pd.Series): predicted energy above convex
            hull
        stability_threshold (float): Where to place stability threshold relative to
            convex hull in eV/atom, usually 0 or 0.1 eV. Default = STABILITY_THRESHOLD.
        fillna (bool): Whether to fill NaNs as the model predicting unstable. Defaults
            to True.

    Note: Should give equivalent classification metrics to
        sklearn.metrics.classification_report(
            each_true > stability_threshold,
            each_pred > stability_threshold,
            output_dict=True,
        )
        when using the same stability_threshold.

    Returns:
        dict[str, float]: dictionary of classification metrics with keys DAF, Precision,
            Recall, Accuracy, F1, TP, FP, TN, FN, MAE, RMSE, R2.
    """
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum,
        classify_stable(
            each_true, each_pred, stability_threshold=stability_threshold, fillna=fillna
        ),
    )

    n_total_pos = n_true_pos + n_false_neg
    n_total_neg = n_true_neg + n_false_pos
    # prevalence: dummy discovery rate of stable crystals by selecting randomly from
    # all materials
    prevalence = (
        n_total_pos / (n_total_pos + n_total_neg)
        if (n_total_pos + n_total_neg) > 0
        else float("nan")
    )
    # Calculate ratios with guards against division by zero
    precision = (
        n_true_pos / (n_true_pos + n_false_pos)
        if (n_true_pos + n_false_pos) > 0
        else float("nan")
    )
    recall = n_true_pos / n_total_pos if n_total_pos > 0 else float("nan")

    # Drop NaNs to calculate regression metrics
    each_true_arr = pd.to_numeric(pd.Series(each_true), errors="coerce")
    each_pred_arr = pd.to_numeric(pd.Series(each_pred), errors="coerce")
    is_nan = each_true_arr.isna() | each_pred_arr.isna()
    each_true = each_true_arr[~is_nan].to_numpy()
    each_pred = each_pred_arr[~is_nan].to_numpy()

    if precision + recall == 0:  # Calculate F1 score, handling division by zero
        f1_score = float("nan")
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    return dict(
        F1=f1_score,
        DAF=precision / prevalence if prevalence > 0 else float("nan"),
        Precision=precision,
        Recall=recall,
        Accuracy=(
            (n_true_pos + n_true_neg) / (n_total_pos + n_total_neg)
            if (n_total_pos + n_total_neg > 0)
            else float("nan")
        ),
        **dict(TP=n_true_pos, FP=n_false_pos, TN=n_true_neg, FN=n_false_neg),
        MAE=np.abs(each_true - each_pred).mean(),
        RMSE=((each_true - each_pred) ** 2).mean() ** 0.5,
        R2=r2_score(each_true, each_pred) if len(each_true) > 1 else float("nan"),
    )


def _align_preds(df_wbm: pd.DataFrame, model_preds: pd.Series) -> pd.Series:
    """Validate prediction IDs against the reference index and coerce to numeric."""
    if unknown_ids := set(model_preds.index) - set(df_wbm.index):
        raise ValueError(
            f"Predictions contain unknown material IDs: {sorted(unknown_ids)}"
        )
    return pd.to_numeric(model_preds.reindex(df_wbm.index), errors="coerce")


def wbm_uniq_proto_prevalence() -> float:
    """Fraction of stable materials among canonical WBM unique prototypes.

    The DAF denominator, computed from the unrounded hull distances so re-evaluated
    models stay comparable to published leaderboard values.
    """
    from matbench_discovery.data import df_wbm

    each_true_uniq = df_wbm.query(MbdKey.uniq_proto)[MbdKey.each_true]
    return float((each_true_uniq <= STABILITY_THRESHOLD).mean())


def discovery_subset_indices(df_wbm: pd.DataFrame) -> dict[TestSubset, pd.Index]:
    """Return canonical WBM subset indices."""
    return {
        TestSubset.full_test_set: df_wbm.index,
        TestSubset.uniq_protos: df_wbm.index[df_wbm[MbdKey.uniq_proto].astype(bool)],
    }


def prepare_model_predictions(
    df_reference: pd.DataFrame,
    model_preds: pd.Series,
    *,
    max_error_threshold: float = MAX_E_FORM_ERROR_THRESHOLD,
) -> tuple[pd.DataFrame, pd.Series]:
    """Clean and align discovery predictions using the leaderboard convention.

    Predictions more than ``max_error_threshold`` eV/atom from DFT are masked.
    Reference columns and predictions are then rounded to three decimals.
    """
    if max_error_threshold < 0:
        raise ValueError(f"{max_error_threshold=} must be nonnegative")
    predictions = pd.to_numeric(
        model_preds.reindex(df_reference.index), errors="coerce"
    )
    bad_mask = abs(predictions - df_reference[MbdKey.e_form_dft]) > max_error_threshold
    predictions = predictions.mask(bad_mask).round(3)
    metric_reference = df_reference.loc[
        :, [MbdKey.each_true, MbdKey.e_form_dft, MbdKey.uniq_proto]
    ].round(3)
    return metric_reference, predictions


def calc_discovery_metrics(
    df_wbm: pd.DataFrame,
    model_preds: pd.Series,
    *,
    uniq_proto_prevalence: float | None = None,
) -> dict[TestSubset, dict[str, float]]:
    """Calculate discovery metrics for both canonical WBM test subsets.

    ``model_preds`` contains formation energies in eV/atom. Predicted hull distances
    use the fixed DFT convex hull, matching the leaderboard and eval script. Reference
    columns and model predictions must use the same rounding convention.

    ``uniq_proto_prevalence`` is the DAF denominator for the uniq-proto subset.
    Callers evaluating against the canonical WBM test set must pass
    :func:`wbm_uniq_proto_prevalence` since ``df_wbm`` reference columns are
    conventionally rounded to 3 decimals, which flips ~430 barely-unstable unique
    prototypes to stable and would silently inflate the prevalence by ~1.3% relative
    to all published DAF values. Defaults to the prevalence of the (possibly rounded)
    ``df_wbm`` frame, intended for synthetic test data only.
    """
    required_cols = {
        str(MbdKey.each_true),
        str(MbdKey.e_form_dft),
        str(MbdKey.uniq_proto),
    }
    if missing_cols := required_cols - set(df_wbm):
        raise ValueError(f"WBM dataframe missing columns: {sorted(missing_cols)}")

    model_preds = _align_preds(df_wbm, model_preds)
    each_true = df_wbm[MbdKey.each_true]
    each_pred = each_true + model_preds - df_wbm[MbdKey.e_form_dft]
    subset_indices = discovery_subset_indices(df_wbm)
    metrics_by_subset = {
        subset: stable_metrics(
            each_true.loc[subset_idx], each_pred.loc[subset_idx], fillna=True
        )
        for subset, subset_idx in subset_indices.items()
    }

    if uniq_proto_prevalence is None:
        each_true_uniq = each_true.loc[subset_indices[TestSubset.uniq_protos]]
        uniq_proto_prevalence = (each_true_uniq <= STABILITY_THRESHOLD).mean()
    daf_denominator = (
        uniq_proto_prevalence if uniq_proto_prevalence > 0 else float("nan")
    )
    uniq_proto_metrics = metrics_by_subset[TestSubset.uniq_protos]
    uniq_proto_metrics["DAF"] = uniq_proto_metrics["Precision"] / daf_denominator
    return metrics_by_subset


def write_all_metrics_to_yaml(
    model: Model,
    metrics_by_subset: Mapping[TestSubset, Mapping[str, float]],
    df_wbm: pd.DataFrame,
    model_preds: pd.Series,
) -> dict[TestSubset, dict[str, str | float]]:
    """Round and rewrite metrics.discovery in one locked update.

    Replaces every subset block (dropping deprecated rate keys) and removes
    obsolete siblings like most_stable_10k; keeps pred_file and cost provenance.
    """
    from ruamel.yaml.comments import CommentedMap

    from matbench_discovery.data import update_yaml_file

    units = {
        "MAE": "eV/atom",
        "RMSE": "eV/atom",
        "R2": "dimensionless",
        "DAF": "dimensionless",
        **dict.fromkeys(("Precision", "Recall", "Accuracy", "F1"), "fraction"),
        **dict.fromkeys(("TP", "FP", "TN", "FN", str(MbdKey.missing_preds)), "count"),
    }
    model_preds = _align_preds(df_wbm, model_preds)
    subset_indices = discovery_subset_indices(df_wbm)
    written: dict[TestSubset, dict[str, str | float]] = {}
    for test_subset, metrics in metrics_by_subset.items():
        block = CommentedMap(
            {key: round(float(value), 3) for key, value in metrics.items()}
        )
        block[str(MbdKey.missing_preds)] = int(
            model_preds.reindex(subset_indices[test_subset]).isna().sum()
        )
        for key, unit in units.items():
            if key in block:
                block.yaml_add_eol_comment(unit, key, column=1)
        written[test_subset] = block

    # Metric subset blocks always carry F1; pred_file / hardware do not.
    update_yaml_file(
        model.yaml_path,
        "metrics.discovery",
        lambda section: {
            **{
                key: val
                for key, val in section.items()
                if not (isinstance(val, Mapping) and "F1" in val)
            },
            **{str(subset): block for subset, block in written.items()},
        },
    )
    return written
