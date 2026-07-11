"""Tests for discovery classification and regression metrics."""

import math
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics import discovery as discovery_module
from matbench_discovery.metrics import metrics_df_from_yaml
from matbench_discovery.metrics.discovery import (
    calc_discovery_metrics,
    classify_stable,
    discovery_subset_indices,
    stable_metrics,
    write_metrics_to_yaml,
)
from matbench_discovery.preds.discovery import df_each_err, df_each_pred, df_preds

df_full_discovery_metrics = metrics_df_from_yaml(["discovery.full_test_set"])
full_discovery_model_labels = set(df_full_discovery_metrics.index)


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

    # Test nullable prediction inputs promised by the type hint
    result = classify_stable([0.0, -0.1], [None, 0.1], stability_threshold=0.0)
    assert [sum(x) for x in result] == [0, 2, 0, 0]


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

    nullable_metrics = stable_metrics(
        [0.0, None, -0.1], [0.0, 0.1, None], stability_threshold=0.0
    )
    assert nullable_metrics["FN"] == 1
    assert nullable_metrics["MAE"] == 0


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


def test_calc_discovery_metrics_matches_manual_three_subset_calculation() -> None:
    """Reusable three-subset metrics preserve the former eval-script semantics."""
    material_ids = pd.Index([f"wbm-{idx}" for idx in range(6)])
    df_test = pd.DataFrame(
        {
            MbdKey.each_true: [-0.2, -0.1, 0.1, 0.2, -0.05, 0.3],
            MbdKey.e_form_dft: [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5],
            MbdKey.uniq_proto: [True, True, False, True, True, False],
        },
        index=material_ids,
    )
    model_preds = pd.Series([-1.1, -0.7, -0.9, -0.6, np.nan, -0.4], index=material_ids)
    metrics_by_subset = calc_discovery_metrics(df_test, model_preds)
    subset_indices = discovery_subset_indices(df_test, model_preds)
    each_pred = df_test[MbdKey.each_true] + model_preds - df_test[MbdKey.e_form_dft]

    assert set(metrics_by_subset) == set(TestSubset)
    assert list(subset_indices[TestSubset.most_stable_10k]) == [
        "wbm-0",
        "wbm-1",
        "wbm-3",
        "wbm-4",
    ]
    assert metrics_by_subset[TestSubset.most_stable_10k]["DAF"] == pytest.approx(4 / 3)
    for subset, subset_idx in subset_indices.items():
        expected = stable_metrics(
            df_test.loc[subset_idx, MbdKey.each_true],
            each_pred.loc[subset_idx],
            fillna=True,
        )
        if subset == TestSubset.full_test_set:
            assert metrics_by_subset[subset] == pytest.approx(expected, nan_ok=True)
        else:
            uniq_prevalence = (
                df_test.loc[subset_indices[TestSubset.uniq_protos], MbdKey.each_true]
                <= 0
            ).mean()
            expected["DAF"] = expected["Precision"] / uniq_prevalence
            assert metrics_by_subset[subset] == pytest.approx(expected, nan_ok=True)


def test_precomputed_discovery_subset_indices_are_reused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Metric calculation and YAML writing reuse one normalized subset ranking."""
    material_ids = pd.Index([f"wbm-{idx}" for idx in range(4)])
    df_test = pd.DataFrame(
        {
            MbdKey.each_true: [-0.2, -0.1, 0.1, 0.2],
            MbdKey.e_form_dft: [-1.0, -0.9, -0.8, -0.7],
            MbdKey.uniq_proto: [True, True, False, True],
        },
        index=material_ids,
    )
    model_preds = pd.to_numeric(
        pd.Series(["-1.1", "-0.7", None, "-0.6"], index=material_ids),
        errors="coerce",
    )
    subset_indices = discovery_subset_indices(df_test, model_preds)

    def reject_recomputation(
        _df_wbm: pd.DataFrame, _model_preds: pd.Series
    ) -> dict[TestSubset, pd.Index]:
        """Fail if either public operation recalculates subset rankings."""
        raise AssertionError("subset rankings were recomputed")

    monkeypatch.setattr(
        discovery_module, "discovery_subset_indices", reject_recomputation
    )
    metrics_by_subset = discovery_module.calc_discovery_metrics(
        df_test, model_preds, subset_indices=subset_indices
    )
    test_yaml = tmp_path / "test_model.yml"
    test_yaml.write_text("metrics:\n  discovery: {}\n")
    mock_model = cast("Model", SimpleNamespace(yaml_path=str(test_yaml)))
    written_metrics = discovery_module.write_all_metrics_to_yaml(
        mock_model,
        metrics_by_subset,
        df_test,
        model_preds,
        subset_indices=subset_indices,
    )
    assert set(written_metrics) == set(TestSubset)
    assert written_metrics[TestSubset.full_test_set][str(MbdKey.missing_preds)] == 1


def test_df_discovery_metrics() -> None:
    """Discovery metrics stored in model YAML files are complete and valid."""
    assert {model.label for model in Model.active()} <= full_discovery_model_labels
    assert df_full_discovery_metrics.MAE.between(0, 0.2).all()
    assert df_full_discovery_metrics.R2.between(-1.5, 1).all()
    assert df_full_discovery_metrics.RMSE.between(0, 0.3).all()
    assert df_full_discovery_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_discovery_eval_skips_incomplete_cli_model() -> None:
    """Test discovery eval skips incomplete CLI models before writing metrics."""
    result = subprocess.run(
        [
            sys.executable,
            "scripts/evals/discovery.py",
            "--models",
            "alphanet-mptrj",
        ],
        cwd=f"{Path(__file__).parents[2]}",
        env=os.environ | {"CI": "1"},
        text=True,
        capture_output=True,
        check=False,
    )
    output = f"{result.stdout}\n{result.stderr}"

    assert result.returncode == 0, output
    assert "Skipping AlphaNet-v1-MPtrj: incomplete discovery metrics" in output
    assert "Loading preds" not in output
    assert "Error processing" not in output
    assert "KeyError" not in output


def test_write_metrics_to_yaml(tmp_path: Path) -> None:
    """Test write_metrics_to_yaml writes metrics with comments to YAML file."""
    test_yaml = tmp_path / "test_model.yml"
    test_yaml.write_text("metrics:\n  discovery: {}\n")
    mock_model = cast("Model", SimpleNamespace(yaml_path=str(test_yaml)))

    # Create test metrics and predictions
    test_metrics: dict[str, str | float] = {
        "MAE": 0.05,
        "RMSE": 0.08,
        "Precision": 0.9,
        "Recall": 0.85,
    }
    test_preds = pd.Series([0.1, 0.2, np.nan, 0.4])  # 1 missing prediction

    result = write_metrics_to_yaml(
        mock_model,
        test_metrics,
        test_preds,
        TestSubset.full_test_set,
    )

    # Check that missing_preds was added
    assert str(MbdKey.missing_preds) in result
    assert result[str(MbdKey.missing_preds)] == 1

    # Check original metrics are preserved
    assert result["MAE"] == 0.05
    assert result["Precision"] == 0.9

    # Verify YAML file was updated
    content = test_yaml.read_text()
    assert "MAE" in content
    assert "full_test_set" in content


def test_df_each_pred() -> None:
    """Assembled per-model predictions span all WBM rows and the metric columns."""
    n_rows, _n_cols = df_each_pred.shape
    assert n_rows == len(df_wbm)
    assert n_rows == len(df_preds)
    assert set(df_each_pred) == full_discovery_model_labels, (
        f"{df_each_pred.columns=}, expected {full_discovery_model_labels=}"
    )


def test_df_each_err() -> None:
    """Per-model hull-distance errors span all WBM rows plus each_err_models."""
    assert len(df_each_err) == len(df_wbm)
    assert len(df_each_err) == len(df_preds)
    assert set(df_each_err) == full_discovery_model_labels | {MbdKey.each_err_models}, (
        f"{df_each_err.columns=}, expected {full_discovery_model_labels=}"
    )
