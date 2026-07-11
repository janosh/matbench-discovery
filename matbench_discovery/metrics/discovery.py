"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.enums import MbdKey, Model, TestSubset

__author__ = "Janosh Riebesell"
__date__ = "2023-02-01"


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
            Recall, Accuracy, F1, TPR, FPR, TNR, FNR, MAE, RMSE, R2.

    Raises:
        ValueError: If FPR + TNR don't add up to 1.
        ValueError: If TPR + FNR don't add up to 1.
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

    TPR = recall  # noqa: N806
    FPR = n_false_pos / n_total_neg if n_total_neg > 0 else float("nan")  # noqa: N806
    TNR = n_true_neg / n_total_neg if n_total_neg > 0 else float("nan")  # noqa: N806
    FNR = n_false_neg / n_total_pos if n_total_pos > 0 else float("nan")  # noqa: N806

    # sanity check: false positives + true negatives = all negatives
    if FPR > 0 and TNR > 0 and FPR + TNR != 1:
        raise ValueError(f"{FPR=} {TNR=} don't add up to 1")

    # sanity check: true positives + false negatives = all positives
    if TPR > 0 and FNR > 0 and TPR + FNR != 1:
        raise ValueError(f"{TPR=} {FNR=} don't add up to 1")

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
        **dict(TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR),
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


def discovery_subset_indices(
    df_wbm: pd.DataFrame, model_preds: pd.Series
) -> dict[TestSubset, pd.Index]:
    """Return canonical WBM subset indices.

    The most-stable 10k are ranked by predicted hull distance, not raw formation
    energy, because discovery uses the fixed DFT convex hull.
    """
    model_preds = _align_preds(df_wbm, model_preds)
    each_pred = df_wbm[MbdKey.each_true] + model_preds - df_wbm[MbdKey.e_form_dft]
    uniq_proto_idx = df_wbm.index[df_wbm[MbdKey.uniq_proto].astype(bool)]
    most_stable_10k_idx = (
        each_pred.loc[uniq_proto_idx]
        .sort_values(na_position="last", kind="stable")
        .head(10_000)
        .index
    )
    return {
        TestSubset.full_test_set: df_wbm.index,
        TestSubset.uniq_protos: uniq_proto_idx,
        TestSubset.most_stable_10k: most_stable_10k_idx,
    }


def calc_discovery_metrics(
    df_wbm: pd.DataFrame,
    model_preds: pd.Series,
    *,
    subset_indices: Mapping[TestSubset, pd.Index] | None = None,
    uniq_proto_prevalence: float | None = None,
) -> dict[TestSubset, dict[str, float]]:
    """Calculate discovery metrics for all three canonical WBM test subsets.

    ``model_preds`` contains formation energies in eV/atom. Predicted hull distances
    use the fixed DFT convex hull, matching the leaderboard and eval script. Reference
    columns and model predictions must use the same rounding convention.
    Pass ``subset_indices`` to reuse rankings already derived from these predictions.

    ``uniq_proto_prevalence`` is the DAF denominator for the uniq-proto and 10k
    subsets. Callers evaluating against the canonical WBM test set must pass
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
    if subset_indices is None:
        subset_indices = discovery_subset_indices(df_wbm, model_preds)
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
    for subset in (TestSubset.uniq_protos, TestSubset.most_stable_10k):
        metrics_by_subset[subset]["DAF"] = (
            metrics_by_subset[subset]["Precision"] / daf_denominator
        )
    return metrics_by_subset


def write_all_metrics_to_yaml(
    model: Model,
    metrics_by_subset: Mapping[TestSubset, Mapping[str, float]],
    df_wbm: pd.DataFrame,
    model_preds: pd.Series,
    *,
    subset_indices: Mapping[TestSubset, pd.Index] | None = None,
) -> dict[TestSubset, dict[str, str | float]]:
    """Round and write all canonical discovery subsets to one model YAML.

    Pass the ``subset_indices`` from :func:`discovery_subset_indices` that also fed
    :func:`calc_discovery_metrics` to avoid reranking predictions before writing
    subset predictions and missing counts.
    """
    model_preds = _align_preds(df_wbm, model_preds)
    if subset_indices is None:
        subset_indices = discovery_subset_indices(df_wbm, model_preds)
    return {
        test_subset: write_metrics_to_yaml(
            model,
            {key: round(float(value), 3) for key, value in metrics.items()},
            model_preds.reindex(subset_indices[test_subset]),
            test_subset,
        )
        for test_subset, metrics in metrics_by_subset.items()
    }


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, str | float],
    df_model_preds: pd.Series,
    test_subset: TestSubset,
) -> dict[str, str | float]:
    """Write discovery metrics to model's YAML file.

    Args:
        model (Model): Model to write metrics for.
        metrics (dict[str, float]): Metrics for this model and test subset.
        df_model_preds (pd.Series): Model predictions for this test subset.
        test_subset (TestSubset): Which test subset these metrics are for.

    Returns:
        dict[str, str | float]: Discovery metrics for this model and test subset.
    """
    from ruamel.yaml.comments import CommentedMap

    from matbench_discovery.data import update_yaml_file

    # calculate number of missing predictions
    n_missing = int(df_model_preds.isna().sum())
    metrics[str(MbdKey.missing_preds)] = n_missing

    # Define metric units for end-of-line comments
    metric_units = {
        "MAE": "eV/atom",
        "RMSE": "eV/atom",
        "R2": "dimensionless",
        "DAF": "dimensionless",
        "Precision": "fraction",
        "Recall": "fraction",
        "Accuracy": "fraction",
        "F1": "fraction",
        "TPR": "fraction",
        "FPR": "fraction",
        "TNR": "fraction",
        "FNR": "fraction",
        str(MbdKey.missing_preds): "count",
        "TP": "count",
        "FP": "count",
        "TN": "count",
        "FN": "count",
    }

    # Create CommentedMap and add units as end-of-line comments
    commented_metrics = CommentedMap(metrics)
    for key in metrics:
        if unit := metric_units.get(key):
            commented_metrics.yaml_add_eol_comment(unit, key, column=1)

    # Write back to file
    update_yaml_file(
        model.yaml_path, f"metrics.discovery.{test_subset}", commented_metrics
    )

    return commented_metrics
