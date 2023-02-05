from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

"""Functions to classify energy above convex hull predictions as true/false
positive/negative and compute performance metrics.
"""

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
        stability_threshold (float | None, optional): Maximum energy above convex hull for a
            material to still be considered stable. Usually 0, 0.05 or 0.1. Defaults to
            0, meaning a material has to be directly on the hull to be called stable.
            Negative values mean a material has to pull the known hull down by that
            amount to count as stable. Few materials lie below the known hull, so only
            negative values very close to 0 make sense.

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
    true: Sequence[float], pred: Sequence[float], stability_threshold: float = 0
) -> dict[str, float]:
    """
    Get a dictionary of stability prediction metrics. Mostly binary classification
    metrics, but also MAE, RMSE and R2.

    Args:
        true (list[float]): true energy values
        pred (list[float]): predicted energy values
        stability_threshold (float): Where to place stability threshold relative to
            convex hull in eV/atom, usually 0 or 0.1 eV. Defaults to 0.

    Note: Could be replaced by sklearn.metrics.classification_report() which takes
        binary labels. I.e. classification_report(true > 0, pred > 0, output_dict=True)
        should give equivalent results.

    Returns:
        dict[str, float]: dictionary of classification metrics with keys DAF, Precision,
            Recall, Accuracy, F1, TPR, FPR, TNR, FNR, MAE, RMSE, R2.
    """
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        true, pred, stability_threshold
    )

    n_true_pos, n_false_pos, n_true_neg, n_false_neg = map(
        sum, (true_pos, false_pos, true_neg, false_neg)
    )

    n_total_pos = n_true_pos + n_false_neg
    n_total_neg = n_true_neg + n_false_pos
    prevalence = n_total_pos / len(true)  # null rate
    precision = n_true_pos / (n_true_pos + n_false_pos)
    recall = n_true_pos / n_total_pos

    is_nan = np.isnan(true) | np.isnan(pred)
    true, pred = np.array(true)[~is_nan], np.array(pred)[~is_nan]

    return dict(
        DAF=precision / prevalence,
        Precision=precision,
        Recall=recall,
        Accuracy=(n_true_pos + n_true_neg) / len(true),
        F1=2 * (precision * recall) / (precision + recall),
        TPR=n_true_pos / n_total_pos,
        FPR=n_false_pos / n_total_neg,
        TNR=n_true_neg / n_total_neg,
        FNR=n_false_neg / n_total_pos,
        MAE=np.abs(true - pred).mean(),
        RMSE=((true - pred) ** 2).mean() ** 0.5,
        R2=r2_score(true, pred),
    )
