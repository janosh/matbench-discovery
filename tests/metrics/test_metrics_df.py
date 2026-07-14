"""Tests for metrics_df_from_yaml aggregation across modeling tasks."""

from matbench_discovery.enums import TestSubset
from matbench_discovery.metrics import metrics_df_from_yaml


def test_metrics_df_from_yaml_discovery() -> None:
    """Discovery unique-prototype metrics load with non-negative MAE."""
    df_metrics = metrics_df_from_yaml([f"discovery.{TestSubset.uniq_protos}"])
    assert df_metrics.shape[0] >= 24
    assert df_metrics.shape[1] >= 15
    assert "MAE" in df_metrics
    assert all(df_metrics.MAE >= 0), f"{df_metrics.MAE=}"
    metrics_df_from_yaml(["discovery.full_test_set"])


def test_metrics_df_from_yaml_phonons() -> None:
    """Phonon kappa metrics load with a single non-negative SRME column."""
    df_metrics = metrics_df_from_yaml(["phonons.kappa_103"])
    kappa_cols = [col for col in df_metrics if col.endswith("_SRME")]
    assert len(kappa_cols) == 1, f"{df_metrics.columns=}"
    assert df_metrics.shape[0] >= 16
    assert df_metrics.shape[1] >= 2
    assert all(df_metrics[kappa_cols[0]] >= 0), f"{df_metrics[kappa_cols[0]]=}"


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
    assert len(df_metrics_both.columns) == len(df_metrics_1.columns) + len(
        df_metrics_2.columns
    )


def test_metrics_df_from_yaml_edge_cases() -> None:
    """Empty or invalid metric paths return an empty frame."""
    assert metrics_df_from_yaml([]).empty
    assert metrics_df_from_yaml(["invalid.path"]).empty
