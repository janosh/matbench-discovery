"""Tests for metrics_df_from_yaml aggregation across modeling tasks."""

import pytest

from matbench_discovery.enums import TestSubset
from matbench_discovery.metrics import metrics_df_from_yaml


@pytest.mark.parametrize(
    ("nested_key", "metric_col"),
    [
        (f"discovery.{TestSubset.uniq_protos}", "MAE"),
        ("phonons.kappa_103", "κ_SRME"),
    ],
)
def test_metrics_df_from_yaml_task_metrics(nested_key: str, metric_col: str) -> None:
    """Task metrics load into a nonempty frame with valid errors."""
    df_metrics = metrics_df_from_yaml([nested_key])
    assert metric_col in df_metrics
    assert df_metrics[metric_col].ge(0).all()


def test_metrics_df_from_yaml_preserves_non_numeric_values() -> None:
    """String and mapping metadata survive metric aggregation."""
    df_metrics = metrics_df_from_yaml(["md"])
    assert "hardware" in df_metrics
    assert "pred_file" in df_metrics
    assert df_metrics.hardware.dropna().map(type).eq(str).all()
    assert df_metrics.pred_file.dropna().map(type).eq(dict).all()


def test_metrics_df_from_yaml_multiple_paths() -> None:
    """Multiple metric paths concatenate columns without dropping either task."""
    df_metrics_1 = metrics_df_from_yaml([f"discovery.{TestSubset.uniq_protos}"])
    df_metrics_2 = metrics_df_from_yaml(["phonons.kappa_103"])
    df_metrics_both = metrics_df_from_yaml(
        [f"discovery.{TestSubset.uniq_protos}", "phonons.kappa_103"]
    )
    assert set(df_metrics_both) == set(df_metrics_1) | set(df_metrics_2)


@pytest.mark.parametrize(
    "nested_keys",
    [(), ("invalid.path",), ("md.run_time_sec",), ("md.run_time_sec.invalid",)],
)
def test_metrics_df_from_yaml_ignores_invalid_paths(
    nested_keys: tuple[str, ...],
) -> None:
    """Missing and scalar metric paths return an empty frame."""
    assert metrics_df_from_yaml(nested_keys).empty
