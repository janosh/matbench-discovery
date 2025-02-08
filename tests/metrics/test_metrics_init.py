"""Test metrics/__init__.py module."""

from matbench_discovery.enums import TestSubset
from matbench_discovery.metrics import metrics_df_from_yaml


def test_metrics_df_from_yaml_discovery() -> None:
    """Test extracting discovery metrics."""
    df_metrics = metrics_df_from_yaml([f"discovery.{TestSubset.uniq_protos}"])
    assert df_metrics.shape >= (24, 18)
    assert "MAE" in df_metrics
    assert all(df_metrics.MAE >= 0), f"{df_metrics.MAE=}"

    metrics_df_from_yaml(["discovery.full_test_set"])


def test_metrics_df_from_yaml_phonons() -> None:
    """Test extracting phonon metrics."""
    df_metrics = metrics_df_from_yaml(["phonons.kappa_103"])
    assert "κ_SRME" in df_metrics
    assert df_metrics.shape >= (16, 3)
    assert all(df_metrics.κ_SRME >= 0), f"{df_metrics.κ_SRME=}"


def test_metrics_df_from_yaml_multiple_paths() -> None:
    """Test extracting metrics from multiple paths."""
    df_metrics_1 = metrics_df_from_yaml([f"discovery.{TestSubset.uniq_protos}"])
    df_metrics_2 = metrics_df_from_yaml(["phonons.kappa_103"])
    df_metrics_both = metrics_df_from_yaml(
        [f"discovery.{TestSubset.uniq_protos}", "phonons.kappa_103"]
    )
    assert df_metrics_both.shape >= (24, 18)
    assert df_metrics_both.shape >= (16, 3)
    n_cols_1 = len(df_metrics_1.columns)
    n_cols_2 = len(df_metrics_2.columns)
    n_cols_both = len(df_metrics_both.columns)
    assert n_cols_both == n_cols_1 + n_cols_2


def test_metrics_df_from_yaml_edge_cases() -> None:
    """Test edge cases: empty keys, invalid paths, not applicable metrics."""
    assert metrics_df_from_yaml([]).empty  # empty keys
    assert metrics_df_from_yaml(["invalid.path"]).empty  # invalid path
