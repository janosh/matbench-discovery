"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

from collections.abc import Sequence

import numpy as np
import pandas as pd
from pymatviz.enums import Key
from sklearn.metrics import r2_score

from matbench_discovery import PDF_FIGS, STABILITY_THRESHOLD
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics import metrics_df_from_yaml

__author__ = "Janosh Riebesell"
__date__ = "2023-02-01"


def classify_stable(
    e_above_hull_true: Sequence[float],
    e_above_hull_pred: Sequence[float],
    *,
    stability_threshold: float | None = 0,
    fillna: bool = True,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Classify model stability predictions as true/false positive/negatives (usually
    w.r.t DFT-ground truth labels). All energies are assumed to be in eV/atom
    (but shouldn't really matter as long as they're consistent).

    Args:
        e_above_hull_true (Sequence[float]): Ground truth energy above hull
            values.
        e_above_hull_pred (Sequence[float]): Model-predicted energy above hull values.
        stability_threshold (float | None, optional): Maximum energy above convex hull
            for a material to still be considered stable. Usually 0, 0.05 or 0.1.
            Defaults to 0, meaning a material has to be directly on the hull to be
            called stable. Negative values mean a material has to pull the known hull
            down by that amount to count as stable. Few materials lie below the known
            hull, so only negative values very close to 0 make sense.
        fillna (bool): Whether to fill NaNs as the model predicting unstable. Defaults
            to True.

    Returns:
        tuple[TP, FN, FP, TN]: Indices as pd.Series for true positives,
            false negatives, false positives and true negatives (in this order).

    Raises:
        ValueError: If sum of positive + negative preds doesn't add up to the total.
    """
    each_true, each_pred = pd.Series(e_above_hull_true), pd.Series(e_above_hull_pred)

    actual_pos = each_true <= (stability_threshold or 0)  # guard against None
    actual_neg = each_true > (stability_threshold or 0)

    model_pos = each_pred <= (stability_threshold or 0)
    model_neg = each_pred > (stability_threshold or 0)

    if fillna:
        nan_mask = np.isnan(e_above_hull_pred)
        # for in both the model's stable and unstable preds, fill NaNs as unstable
        model_pos[nan_mask] = False
        model_neg[nan_mask] = True

        n_pos, n_neg, total = model_pos.sum(), model_neg.sum(), len(e_above_hull_pred)
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
    each_true: Sequence[float],
    each_pred: Sequence[float],
    *,
    stability_threshold: float = STABILITY_THRESHOLD,
    fillna: bool = True,
) -> dict[str, float]:
    """Get a dictionary of stability prediction metrics. Mostly binary classification
    metrics, but also MAE, RMSE and R2.

    Args:
        each_true (list[float]): true energy above convex hull
        each_pred (list[float]): predicted energy above convex hull
        stability_threshold (float): Where to place stability threshold relative to
            convex hull in eV/atom, usually 0 or 0.1 eV. Defaults to 0.
        fillna (bool): Whether to fill NaNs as the model predicting unstable. Defaults
            to True.

    Note: Should give equivalent classification metrics to
        sklearn.metrics.classification_report(
            each_true > STABILITY_THRESHOLD,
            each_pred > STABILITY_THRESHOLD,
            output_dict=True,
        )

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
    prevalence = n_total_pos / (n_total_pos + n_total_neg)
    precision = n_true_pos / (n_true_pos + n_false_pos)  # model's discovery rate
    recall = n_true_pos / n_total_pos

    TPR = recall
    FPR = n_false_pos / n_total_neg
    TNR = n_true_neg / n_total_neg
    FNR = n_false_neg / n_total_pos

    if FPR + TNR != 1:  # sanity check: false positives + true negatives = all negatives
        raise ValueError(f"{FPR=} {TNR=} don't add up to 1")

    if TPR + FNR != 1:  # sanity check: true positives + false negatives = all positives
        raise ValueError(f"{TPR=} {FNR=} don't add up to 1")

    # Drop NaNs to calculate regression metrics
    is_nan = np.isnan(each_true) | np.isnan(each_pred)
    each_true, each_pred = np.array(each_true)[~is_nan], np.array(each_pred)[~is_nan]

    if precision + recall == 0:  # Calculate F1 score, handling division by zero
        f1_score = float("nan")
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    return dict(
        F1=f1_score,
        DAF=precision / prevalence,
        Precision=precision,
        Recall=recall,
        Accuracy=(n_true_pos + n_true_neg) / len(each_true),
        **dict(TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR),
        **dict(TP=n_true_pos, FP=n_false_pos, TN=n_true_neg, FN=n_false_neg),
        MAE=np.abs(each_true - each_pred).mean(),
        RMSE=((each_true - each_pred) ** 2).mean() ** 0.5,
        R2=r2_score(each_true, each_pred),
    )


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


# Create DataFrames with models as rows
df_metrics = (
    metrics_df_from_yaml(["discovery.full_test_set"])
    .sort_values(by=Key.f1.upper(), ascending=False)
    .T
)
df_metrics_10k = (
    metrics_df_from_yaml(["discovery.most_stable_10k"])
    .sort_values(by=Key.f1.upper(), ascending=False)
    .T
)
df_metrics_uniq_protos = (
    metrics_df_from_yaml(["discovery.unique_prototypes", "phonons"])
    .sort_values(by=Key.f1.upper(), ascending=False)
    .T
)
df_metrics_uniq_protos = df_metrics_uniq_protos.drop(index=[MbdKey.missing_preds])

for df, title in (
    (df_metrics, "Metrics for Full Test Set"),
    (df_metrics_10k, "Metrics for 10k Most Stable Predictions"),
    (df_metrics_uniq_protos, "Metrics for unique non-MP prototypes"),
):
    df.attrs["title"] = title


dfs_metrics: dict[TestSubset, pd.DataFrame] = {
    TestSubset.full_test_set: df_metrics,
    TestSubset.uniq_protos: df_metrics_uniq_protos,
    TestSubset.most_stable_10k: df_metrics_10k,
}

if __name__ == "__main__":
    # export metrics tables to latex
    for test_subset, df_subset in dfs_metrics.items():
        tex_path = f"{PDF_FIGS}/metrics-table-{test_subset.replace('_', '-')}.tex"
        df_subset.to_latex(
            tex_path,
            float_format=lambda x: f"{x:.3f}",
            escape=True,
            bold_rows=True,
            caption=df_subset.attrs["title"],
        )
        print(f"Wrote {tex_path}")
