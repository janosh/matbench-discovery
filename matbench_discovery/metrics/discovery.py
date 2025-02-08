"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

from collections.abc import Sequence

import numpy as np
import pandas as pd
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap
from sklearn.metrics import r2_score

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.data import df_wbm, round_trip_yaml
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics import metrics_df_from_yaml

__author__ = "Janosh Riebesell"
__date__ = "2023-02-01"


def classify_stable(
    e_above_hull_true: pd.Series,
    e_above_hull_pred: pd.Series,
    *,
    stability_threshold: float | None = 0,
    fillna: bool = True,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Classify model stability predictions as true/false positive/negatives (usually
    w.r.t DFT-ground truth labels). All energies are assumed to be in eV/atom
    (but shouldn't really matter as long as they're consistent).

    Args:
        e_above_hull_true (pd.Series): Ground truth energy above convex hull values.
        e_above_hull_pred (pd.Series): Model predicted energy above convex hull values.
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
    actual_pos = e_above_hull_true <= (stability_threshold or 0)  # guard against None
    actual_neg = e_above_hull_true > (stability_threshold or 0)

    model_pos = e_above_hull_pred <= (stability_threshold or 0)
    model_neg = e_above_hull_pred > (stability_threshold or 0)

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

    return dict(
        F1=2 * (precision * recall) / (precision + recall),
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


def write_discovery_metrics_to_yaml(
    model: Model,
    df_metrics: pd.DataFrame,
    df_metrics_10k: pd.DataFrame,
    df_metrics_uniq_protos: pd.DataFrame,
    df_preds: pd.DataFrame,
) -> None:
    """Write materials discovery metrics to model YAML metadata files."""
    metrics_full_test_set = df_metrics[model.label].to_dict()
    metrics_10k_most_stable = df_metrics_10k[model.label].to_dict()
    metrics_unique_protos = df_metrics_uniq_protos[model.label].to_dict()

    df_uniq_proto_preds = df_preds[df_wbm[MbdKey.uniq_proto]]

    each_pred_uniq_proto = (
        df_uniq_proto_preds[MbdKey.each_true]
        + df_uniq_proto_preds[model.label]
        - df_uniq_proto_preds[MbdKey.e_form_dft]
    )
    most_stable_10k_idx = each_pred_uniq_proto.nsmallest(10_000).index

    # calculate number of missing predictions for each test subset
    for metrics, df_tmp in (
        (metrics_full_test_set, df_preds),
        (metrics_10k_most_stable, df_preds.loc[most_stable_10k_idx]),
        (metrics_unique_protos, df_preds.query(MbdKey.uniq_proto)),
    ):
        metrics[str(MbdKey.missing_preds)] = int(df_tmp[model.label].isna().sum())
        metrics[str(MbdKey.missing_percent)] = (
            f"{metrics[str(MbdKey.missing_preds)] / len(df_tmp):.2%}"
        )

    new_metrics = {
        str(key): CommentedMap(metrics)
        for key, metrics in (
            (TestSubset.full_test_set, metrics_full_test_set),
            (TestSubset.most_stable_10k, metrics_10k_most_stable),
            (TestSubset.uniq_protos, metrics_unique_protos),
        )
    }

    # Add or update discovery metrics
    with open(model.yaml_path) as file:
        model_metadata = round_trip_yaml.load(file)

    all_metrics = model_metadata.setdefault("metrics", {})
    discovery_metrics = all_metrics.setdefault("discovery", {})
    discovery_metrics |= new_metrics

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
        str(MbdKey.missing_percent): "fraction",
        str(MbdKey.missing_preds): "count",
        "TP": "count",
        "FP": "count",
        "TN": "count",
        "FN": "count",
    }

    # Add units as YAML end-of-line comments for each test subset
    for key in new_metrics:
        subset_dict: CommentedMap = discovery_metrics[key]
        for key in subset_dict:
            if unit := metric_units.get(key):
                subset_dict.yaml_add_eol_comment(unit, key)

    all_metrics["discovery"] = discovery_metrics

    # Write back to file
    with open(model.yaml_path, mode="w") as file:
        round_trip_yaml.dump(model_metadata, file)


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
df_metrics_uniq_protos = df_metrics_uniq_protos.drop(
    index=[MbdKey.missing_preds, MbdKey.missing_percent]
)

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
