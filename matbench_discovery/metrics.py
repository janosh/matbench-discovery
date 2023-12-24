"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sklearn.metrics import r2_score

from matbench_discovery import STABILITY_THRESHOLD

if TYPE_CHECKING:
    from collections.abc import Sequence

    import pandas as pd


__author__ = "Janosh Riebesell"
__date__ = "2023-02-01"


def classify_stable(
    e_above_hull_true: pd.Series,
    e_above_hull_pred: pd.Series,
    stability_threshold: float | None = 0,
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

    Returns:
        tuple[TP, FN, FP, TN]: Indices as pd.Series for true positives,
            false negatives, false positives and true negatives (in this order).
    """
    actual_pos = e_above_hull_true <= (stability_threshold or 0)  # guard against None
    actual_neg = e_above_hull_true > (stability_threshold or 0)
    model_pos = e_above_hull_pred <= (stability_threshold or 0)
    model_neg = e_above_hull_pred > (stability_threshold or 0)

    true_pos = actual_pos & model_pos
    false_neg = actual_pos & model_neg
    false_pos = actual_neg & model_pos
    true_neg = actual_neg & model_neg

    return true_pos, false_neg, false_pos, true_neg


def stable_metrics(
    each_true: Sequence[float],
    each_pred: Sequence[float],
    stability_threshold: float = STABILITY_THRESHOLD,
) -> dict[str, float]:
    """Get a dictionary of stability prediction metrics. Mostly binary classification
    metrics, but also MAE, RMSE and R2.

    Args:
        each_true (list[float]): true energy above convex hull
        each_pred (list[float]): predicted energy above convex hull
        stability_threshold (float): Where to place stability threshold relative to
            convex hull in eV/atom, usually 0 or 0.1 eV. Defaults to 0.

    Note: Should give equivalent classification metrics to
        sklearn.metrics.classification_report(
            each_true > STABILITY_THRESHOLD,
            each_pred > STABILITY_THRESHOLD,
            output_dict=True,
        )

    Returns:
        dict[str, float]: dictionary of classification metrics with keys DAF, Precision,
            Recall, Accuracy, F1, TPR, FPR, TNR, FNR, MAE, RMSE, R2.
    """
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, classify_stable(each_true, each_pred, stability_threshold)
    )

    n_total_pos = n_true_pos + n_false_neg
    n_total_neg = n_true_neg + n_false_pos
    # prevalence: dummy discovery rate of stable crystals by selecting randomly from
    # all materials
    prevalence = n_total_pos / len(each_true)
    precision = n_true_pos / (n_true_pos + n_false_pos)  # model's discovery rate
    recall = n_true_pos / n_total_pos

    is_nan = np.isnan(each_true) | np.isnan(each_pred)
    each_true, each_pred = np.array(each_true)[~is_nan], np.array(each_pred)[~is_nan]

    TPR = recall
    FPR = n_false_pos / n_total_neg
    TNR = n_true_neg / n_total_neg
    FNR = n_false_neg / n_total_pos

    if FPR + TNR != 1:  # sanity check: false positives + true negatives = all negatives
        raise ValueError(f"{FPR=} {TNR=} don't add up to 1")

    if TPR + FNR != 1:  # sanity check: true positives + false negatives = all positives
        raise ValueError(f"{TPR=} {FNR=} don't add up to 1")

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
