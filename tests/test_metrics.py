from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from pytest import approx

from matbench_discovery.metrics import classify_stable, stable_metrics


@pytest.mark.parametrize(
    "stability_threshold, expected",
    [(-0.1, (6, 10, 7, 7)), (0, (6, 10, 7, 7)), (0.1, (10, 8, 6, 6))],
)
def test_classify_stable(
    stability_threshold: float, expected: tuple[int, int, int, int]
) -> None:
    df = pd._testing.makeDataFrame()

    true_pos, false_neg, false_pos, true_neg = classify_stable(
        e_above_hull_true=df.A,
        e_above_hull_pred=df.B,
        stability_threshold=stability_threshold,
    )
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, (true_pos, false_neg, false_pos, true_neg)
    )

    assert (n_true_pos, n_false_neg, n_false_pos, n_true_neg) == expected
    assert n_true_pos + n_false_neg + n_false_pos + n_true_neg == len(df)
    assert n_true_neg + n_false_pos == sum(stability_threshold < df.A)
    assert n_true_pos + n_false_neg == sum(stability_threshold >= df.A)


def test_stable_metrics() -> None:
    metrics = stable_metrics(np.arange(-1, 1, 0.1), np.arange(1, -1, -0.1))
    for key, val in dict(
        DAF=0,
        Precision=0,
        Recall=0,
        Accuracy=0,
        TPR=0,
        FPR=1,
        TNR=0,
        FNR=1,
        MAE=0.999,
        RMSE=1.157,
        R2=-3.030,
    ).items():
        assert metrics[key] == approx(val, abs=1e-3), f"{key=}"

    assert math.isnan(metrics["F1"])

    # test stable_metrics gives the same result as sklearn.metrics.classification_report
    # for random numpy data
    # np.random.seed(0)
    y_true, y_pred = np.random.randn(100, 2).T
    metrics = stable_metrics(y_true, y_pred)

    from sklearn.metrics import classification_report

    skl_report = classification_report(y_true > 0, y_pred > 0, output_dict=True)

    assert metrics["Precision"] == skl_report["False"]["precision"]
    assert metrics["Recall"] == skl_report["False"]["recall"]
    assert metrics["F1"] == skl_report["False"]["f1-score"]
    assert metrics["Accuracy"] == skl_report["accuracy"]

    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    # test regression metrics
    assert metrics["MAE"] == mean_absolute_error(y_true, y_pred)
    assert metrics["RMSE"] == mean_squared_error(y_true, y_pred) ** 0.5
    assert metrics["R2"] == r2_score(y_true, y_pred)

    # test stable_metrics docstring is up to date, all returned metrics should be listed
    assert stable_metrics.__doc__  # for mypy
    assert all(key in stable_metrics.__doc__ for key in metrics)

    # test discovery acceleration factor (DAF)
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, classify_stable(y_true, y_pred)
    )

    dummy_hit_rate = (n_true_pos + n_false_neg) / (
        n_true_pos + n_false_pos + n_false_neg + n_true_neg
    )
    precision = n_true_pos / (n_true_pos + n_false_pos)
    assert metrics["DAF"] == precision / dummy_hit_rate
