"""Tests for discovery classification and regression metrics."""

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import (
    classification_report,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)

from matbench_discovery import ROOT
from matbench_discovery.data import df_wbm, load_discovery_predictions
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics import discovery as discovery_module
from matbench_discovery.metrics import metrics_df_from_yaml
from matbench_discovery.metrics.discovery import (
    calc_discovery_metrics,
    classify_stable,
    discovery_subset_indices,
    stable_metrics,
)

df_preds, df_each_pred, df_each_err = load_discovery_predictions()
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
    """Classification returns aligned boolean masks with expected counts."""
    classifications = classify_stable(
        each_true=df_float.A,
        each_pred=df_float.B,
        stability_threshold=stability_threshold,
    )
    assert all(
        classification.index.equals(df_float.index)
        for classification in classifications
    )
    assert all(classification.dtype == bool for classification in classifications)
    assert tuple(map(sum, classifications)) == expected


def test_classify_stable_rejects_nan_threshold() -> None:
    """A NaN stability threshold is invalid."""
    with pytest.raises(ValueError, match="stability_threshold must be a real number"):
        classify_stable([-0.1, 0.0, 0.1], [-0.1, 0.0, 0.1], stability_threshold=np.nan)


@pytest.mark.parametrize(
    ("each_true", "each_pred", "fillna", "expected"),
    [
        (
            [-0.1, 0.0, 0.1, np.nan, None],
            [-0.1, 0.0, 0.1, 0.2, -0.2],
            True,
            (2, 0, 0, 1),
        ),
        (
            np.array([-0.1, 0.0, 0.1, np.nan]),
            np.array([-0.1, 0.0, 0.1, 0.2]),
            False,
            (2, 0, 0, 1),
        ),
        ([0.0, -0.1], [None, 0.1], True, (0, 2, 0, 0)),
    ],
    ids=["nullable-truth", "numpy-no-fill", "nullable-predictions"],
)
def test_classify_stable_nullable_inputs(
    each_true: list[float | None] | np.ndarray,
    each_pred: list[float | None] | np.ndarray,
    fillna: bool,
    expected: tuple[int, int, int, int],
) -> None:
    """Classification supports nullable list and NumPy inputs."""
    result = classify_stable(
        each_true, each_pred, stability_threshold=0.0, fillna=fillna
    )
    assert tuple(map(sum, result)) == expected


@pytest.mark.parametrize(
    ("each_true", "each_pred", "expected"),
    [
        (
            [0.1, 0.2, 0.3],
            [-0.1, -0.2, -0.3],
            (0.0, 1.0, 0.0, np.nan, np.nan, np.nan, np.nan),
        ),
        (
            [-0.1, -0.2, -0.3],
            [0.1, 0.2, 0.3],
            (np.nan, np.nan, np.nan, 0.0, 1.0, np.nan, np.nan),
        ),
    ],
    ids=["no-actual-stable", "no-actual-unstable"],
)
def test_stable_metrics_zero_class(
    each_true: list[float],
    each_pred: list[float],
    expected: tuple[float, ...],
) -> None:
    """Metrics handle absent stable or unstable classes."""
    metrics = stable_metrics(each_true, each_pred, stability_threshold=0.0, fillna=True)
    metric_keys = ("Precision", "FPR", "TNR", "Recall", "FNR", "DAF", "F1")
    assert tuple(metrics[key] for key in metric_keys) == pytest.approx(
        expected, nan_ok=True
    )


def test_stable_metrics_regression_edge_cases() -> None:
    """Regression metrics handle one or no valid pairs."""
    assert np.isnan(stable_metrics([0.1], [0.2], stability_threshold=0.0)["R2"])
    all_nan_metrics = stable_metrics(
        [np.nan, np.nan], [np.nan, np.nan], stability_threshold=0.0
    )
    assert all(np.isnan(all_nan_metrics[key]) for key in ["MAE", "RMSE", "R2"])


def test_stable_metrics_nan_handling() -> None:
    """NaN filling only changes classification, not regression metrics."""
    true_values = np.array([-1, 1, 0.1, -0.5, 0.5])
    pred_values = np.array([-1, 1, -0.1, np.nan, np.nan])
    metrics_no_fill = stable_metrics(true_values, pred_values, fillna=False)
    metrics_fill = stable_metrics(true_values, pred_values, fillna=True)

    assert metrics_no_fill["Precision"] == metrics_fill["Precision"]
    assert metrics_no_fill["DAF"] > metrics_fill["DAF"]
    assert metrics_no_fill["TNR"] == 0.5
    assert metrics_no_fill["FNR"] == 0
    assert metrics_fill["TNR"] == 2 / 3
    assert metrics_fill["FNR"] == 1 / 2
    for metric in ("MAE", "RMSE", "R2"):
        assert metrics_no_fill[metric] == metrics_fill[metric]

    nullable_metrics = stable_metrics(
        [0.0, None, -0.1], [0.0, 0.1, None], stability_threshold=0.0
    )
    assert nullable_metrics["FN"] == 1
    assert nullable_metrics["MAE"] == 0


def test_stable_metrics_matches_sklearn() -> None:
    """Classification and regression metrics match sklearn."""
    np_rng = np.random.default_rng(seed=0)
    true_values, pred_values = np_rng.normal(size=(2, 100))
    metrics = stable_metrics(true_values, pred_values, fillna=True)
    sklearn_report = classification_report(
        true_values > 0, pred_values > 0, output_dict=True
    )
    expected = {
        "Precision": sklearn_report["False"]["precision"],
        "Recall": sklearn_report["False"]["recall"],
        "F1": sklearn_report["False"]["f1-score"],
        "Accuracy": sklearn_report["accuracy"],
        "MAE": mean_absolute_error(true_values, pred_values),
        "RMSE": mean_squared_error(true_values, pred_values) ** 0.5,
        "R2": r2_score(true_values, pred_values),
    }
    assert {key: metrics[key] for key in expected} == pytest.approx(expected)


def test_most_stable_10k_includes_missing_predictions_last() -> None:
    """Missing predictions fill remaining 10k slots after ranked predictions."""
    material_ids = pd.Index([f"wbm-{idx}" for idx in range(10_001)])
    df_test = pd.DataFrame(
        {MbdKey.each_true: 0.0, MbdKey.e_form_dft: 0.0, MbdKey.uniq_proto: True},
        index=material_ids,
    )
    model_preds = pd.Series(
        [*range(9_999), np.nan, np.nan], index=material_ids, dtype=float
    )

    subset_idx = discovery_subset_indices(df_test, model_preds)[
        TestSubset.most_stable_10k
    ]
    assert len(subset_idx) == 10_000
    assert subset_idx[:-1].equals(material_ids[:9_999])
    assert pd.isna(model_preds.loc[subset_idx[-1]])


def test_most_stable_10k_preserves_unique_prototype_order_for_ties() -> None:
    """Equal predictions preserve the original unique-prototype order."""
    material_ids = pd.Index([f"wbm-{idx}" for idx in range(10_001)])
    df_test = pd.DataFrame(
        {MbdKey.each_true: 0.0, MbdKey.e_form_dft: 0.0, MbdKey.uniq_proto: True},
        index=material_ids,
    )
    model_preds = pd.Series(0.0, index=material_ids)

    subset_idx = discovery_subset_indices(df_test, model_preds)[
        TestSubset.most_stable_10k
    ]
    assert subset_idx.equals(material_ids[:10_000])


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
    wbm_ids = list(subset_indices[TestSubset.most_stable_10k])
    assert wbm_ids == ["wbm-0", "wbm-1", "wbm-3", "wbm-4"]
    assert metrics_by_subset[TestSubset.most_stable_10k]["DAF"] == pytest.approx(4 / 3)
    uniq_prevalence = (
        df_test.loc[subset_indices[TestSubset.uniq_protos], MbdKey.each_true] <= 0
    ).mean()
    for subset, subset_idx in subset_indices.items():
        expected = stable_metrics(
            df_test.loc[subset_idx, MbdKey.each_true],
            each_pred.loc[subset_idx],
            fillna=True,
        )
        if subset != TestSubset.full_test_set:
            expected["DAF"] = expected["Precision"] / uniq_prevalence
        assert metrics_by_subset[subset] == pytest.approx(expected, nan_ok=True)

    # an explicit prevalence (e.g. from unrounded WBM hull distances) overrides the
    # denominator derived from the possibly-rounded reference frame
    overridden = calc_discovery_metrics(df_test, model_preds, uniq_proto_prevalence=0.5)
    for subset in (TestSubset.uniq_protos, TestSubset.most_stable_10k):
        assert overridden[subset]["DAF"] == pytest.approx(
            overridden[subset]["Precision"] / 0.5
        )


def test_precomputed_discovery_subset_indices_are_reused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Metric calculation and YAML writing reuse one normalized subset ranking."""
    material_ids = pd.Index([f"wbm-{idx}" for idx in range(5)])
    df_test = pd.DataFrame(
        {
            MbdKey.each_true: [-0.2, 0.1, -0.1, 0.2, 0.3],
            MbdKey.e_form_dft: [0.0] * 5,
            MbdKey.uniq_proto: [True] * 5,
        },
        index=material_ids,
    )
    model_preds = pd.Series(
        [0.0, 0.0, None, 0.0, None], index=material_ids, dtype=float
    )
    expected_by_subset = {
        TestSubset.full_test_set: (material_ids, (1, 1, 3), 2),
        TestSubset.uniq_protos: (material_ids[[0, 1, 2]], (1, 1, 1), 1),
        TestSubset.most_stable_10k: (material_ids[[0, 1]], (1, 0, 1), 0),
    }
    subset_indices = {
        test_subset: expected[0] for test_subset, expected in expected_by_subset.items()
    }

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
    for test_subset, (_, expected_counts, missing_count) in expected_by_subset.items():
        assert (
            tuple(metrics_by_subset[test_subset][key] for key in ("TP", "FN", "TN"))
            == expected_counts
        )
        assert written_metrics[test_subset][str(MbdKey.missing_preds)] == missing_count
    yaml_content = test_yaml.read_text()
    assert "full_test_set:" in yaml_content
    assert f"{MbdKey.missing_preds}:" in yaml_content


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
        [sys.executable, "scripts/evals/discovery.py", "--models", "alphanet-v1-mptrj"],
        cwd=ROOT,
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


@pytest.mark.parametrize(
    ("df_values", "expected_columns"),
    [
        (df_each_pred, full_discovery_model_labels),
        (df_each_err, full_discovery_model_labels | {str(MbdKey.each_err_models)}),
    ],
    ids=["predictions", "hull-errors"],
)
def test_discovery_prediction_frames(
    df_values: pd.DataFrame, expected_columns: set[str]
) -> None:
    """Prediction frames span WBM and contain the expected model columns."""
    assert len(df_values) == len(df_wbm) == len(df_preds)
    assert set(df_values) == expected_columns
