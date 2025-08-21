import math

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.enums import Model
from matbench_discovery.metrics import discovery
from matbench_discovery.metrics.discovery import classify_stable, stable_metrics


@pytest.mark.parametrize(
    "stability_threshold, expected",
    [(-0.1, (2, 4, 3, 1)), (0, (2, 4, 3, 1)), (0.1, (4, 3, 2, 1))],
)
def test_classify_stable(
    stability_threshold: float,
    expected: tuple[int, int, int, int],
    df_float: pd.DataFrame,
) -> None:
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        each_true=df_float.A,
        each_pred=df_float.B,
        stability_threshold=stability_threshold,
    )
    assert all(
        out.index.equals(df_float.index)
        for out in (true_pos, false_neg, false_pos, true_neg)
    )
    assert all(out.dtype == bool for out in (true_pos, false_neg, false_pos, true_neg))

    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, (true_pos, false_neg, false_pos, true_neg)
    )

    assert (n_true_pos, n_false_neg, n_false_pos, n_true_neg) == expected
    assert n_true_pos + n_false_neg + n_false_pos + n_true_neg == len(df_float)
    assert n_true_neg + n_false_pos == np.sum(stability_threshold < df_float.A)
    assert n_true_pos + n_false_neg == np.sum(stability_threshold >= df_float.A)


def test_classify_stable_edge_cases() -> None:
    """Test edge cases for classify_stable function."""
    # Test with None threshold (defaults to 0.0)
    result = classify_stable([-0.1, 0.0, 0.1], [-0.1, 0.0, 0.1])
    # true_pos, false_neg, false_pos, true_neg
    assert [sum(x) for x in result] == [2, 0, 0, 1]

    # Test with NaN threshold (should raise ValueError)
    with pytest.raises(ValueError, match="stability_threshold must be a real number"):
        classify_stable([-0.1, 0.0, 0.1], [-0.1, 0.0, 0.1], stability_threshold=np.nan)

    # Test with numeric threshold
    result = classify_stable(
        [-0.1, 0.0, 0.1], [-0.1, 0.0, 0.1], stability_threshold=0.05
    )
    assert [sum(x) for x in result] == [2, 0, 0, 1]


def test_classify_stable_input_types() -> None:
    """Test classify_stable with different input types including NaN/None values."""
    # Test with Python lists containing NaN and None
    result = classify_stable(
        [-0.1, 0.0, 0.1, np.nan, None],
        [-0.1, 0.0, 0.1, 0.2, -0.2],
        stability_threshold=0.0,
        fillna=True,
    )
    # With fillna=True, NaN/None treated as unstable
    assert [sum(x) for x in result] == [2, 0, 0, 1]

    # Test with NumPy arrays containing NaN
    result = classify_stable(
        np.array([-0.1, 0.0, 0.1, np.nan]),
        np.array([-0.1, 0.0, 0.1, 0.2]),
        stability_threshold=0.0,
        fillna=False,
    )
    assert [sum(x) for x in result] == [2, 0, 0, 1]  # With fillna=False, NaN preserved


def test_stable_metrics_edge_cases() -> None:
    """Test edge cases for stable_metrics function."""
    # Test with all negative predictions (zero positives)
    metrics = stable_metrics(
        [0.1, 0.2, 0.3], [-0.1, -0.2, -0.3], stability_threshold=0.0
    )
    assert metrics["Precision"] == 0.0
    assert metrics["FPR"] == 1.0
    assert metrics["TNR"] == 0.0
    assert np.isnan(metrics["Recall"])
    assert np.isnan(metrics["FNR"])

    assert np.isnan(metrics["DAF"])

    # Test with all positive predictions (zero negatives)
    metrics = stable_metrics(
        [-0.1, -0.2, -0.3], [0.1, 0.2, 0.3], stability_threshold=0.0
    )
    assert np.isnan(metrics["Precision"])
    assert np.isnan(metrics["FPR"])
    assert np.isnan(metrics["TNR"])
    assert metrics["Recall"] == 0.0
    assert metrics["FNR"] == 1.0
    assert np.isnan(metrics["DAF"])

    # Test with single data point and all NaN inputs
    assert np.isnan(stable_metrics([0.1], [0.2], stability_threshold=0.0)["R2"])
    all_nan_metrics = stable_metrics(
        [np.nan, np.nan], [np.nan, np.nan], stability_threshold=0.0
    )
    assert all(np.isnan(all_nan_metrics[key]) for key in ["MAE", "RMSE", "R2"])


def test_stable_metrics_nan_handling() -> None:
    """Test stable_metrics with various NaN handling scenarios."""
    true_vals, pred_vals = [0.1, -0.1, 0.2, -0.2], [0.1, -0.1, np.nan, np.nan]

    metrics_fillna = stable_metrics(
        true_vals, pred_vals, stability_threshold=0.0, fillna=True
    )
    metrics_no_fillna = stable_metrics(
        true_vals, pred_vals, stability_threshold=0.0, fillna=False
    )

    # Classification metrics differ due to NaN handling
    assert metrics_fillna["Recall"] != metrics_no_fillna["Recall"]

    # Regression metrics same (NaN values dropped in both cases)
    assert metrics_fillna["MAE"] == metrics_no_fillna["MAE"]
    assert metrics_fillna["RMSE"] == metrics_no_fillna["RMSE"]


def test_stable_metrics() -> None:
    metrics = stable_metrics(np.arange(-1, 1, 0.1), np.arange(1, -1, -0.1), fillna=True)
    expected = dict(
        DAF=0.0, Precision=0, Recall=0, Accuracy=0, TPR=0, FPR=1, TNR=0, FNR=1
    )
    expected |= dict(MAE=0.999, RMSE=1.157, R2=-3.030)
    assert {k: metrics[k] for k in expected} == pytest.approx(expected, abs=1e-3)

    assert math.isnan(metrics["F1"])

    metrics = stable_metrics(
        np.array((-1, 1, 0.1, -0.5, 0.5)),
        np.array((-1, 1, -0.1, np.nan, np.nan)),
        fillna=False,
    )
    fillna_metrics = stable_metrics(
        np.array((-1, 1, 0.1, -0.5, 0.5)),
        np.array((-1, 1, -0.1, np.nan, np.nan)),
        fillna=True,
    )

    assert metrics["Precision"] == fillna_metrics["Precision"]
    assert metrics["DAF"] > fillna_metrics["DAF"]  # nan's dropped in prevalence
    assert metrics["TNR"] == 0.5
    assert metrics["FNR"] == 0
    assert fillna_metrics["TNR"] == 2 / 3
    assert fillna_metrics["FNR"] == 1 / 2

    # test stable_metrics gives the same result as sklearn.metrics.classification_report
    # for random numpy data
    rng = np.random.default_rng(seed=0)
    y_true, y_pred = rng.normal(size=(2, 100))
    metrics = stable_metrics(y_true, y_pred, fillna=True)

    from sklearn.metrics import classification_report

    skl_report = classification_report(y_true > 0, y_pred > 0, output_dict=True)

    assert metrics["Precision"] == pytest.approx(skl_report["False"]["precision"])
    assert metrics["Recall"] == pytest.approx(skl_report["False"]["recall"])
    assert metrics["F1"] == pytest.approx(skl_report["False"]["f1-score"])
    assert metrics["Accuracy"] == pytest.approx(skl_report["accuracy"])

    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    # test regression metrics
    assert metrics["MAE"] == mean_absolute_error(y_true, y_pred)
    assert metrics["RMSE"] == mean_squared_error(y_true, y_pred) ** 0.5
    assert metrics["R2"] == r2_score(y_true, y_pred)

    # test stable_metrics docstring is up to date, all returned metrics should be listed
    assert stable_metrics.__doc__  # for ty
    assert all(key in stable_metrics.__doc__ for key in metrics)

    # test discovery acceleration factor (DAF)
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, classify_stable(y_true, y_pred, fillna=True)
    )

    dummy_hit_rate = (n_true_pos + n_false_neg) / (
        n_true_pos + n_false_pos + n_false_neg + n_true_neg
    )
    precision = n_true_pos / (n_true_pos + n_false_pos)
    assert metrics[str(Key.daf.symbol)] == precision / dummy_hit_rate


def test_df_discovery_metrics() -> None:
    missing_cols = {*discovery.df_metrics} - {model.label for model in Model}
    assert missing_cols == set(), f"{missing_cols=}"
    assert discovery.df_metrics.T.MAE.between(0, 0.2).all(), (
        f"unexpected {discovery.df_metrics.T.MAE=}"
    )
    assert discovery.df_metrics.T.R2.between(-1.5, 1).all(), (
        f"unexpected {discovery.df_metrics.T.R2=}"
    )
    assert discovery.df_metrics.T.RMSE.between(0, 0.3).all(), (
        f"unexpected {discovery.df_metrics.T.RMSE=}"
    )
    assert discovery.df_metrics.T.isna().sum().sum() == 0, "NaNs in metrics"
