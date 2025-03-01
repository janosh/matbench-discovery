from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.discovery import df_metrics
from matbench_discovery.preds import discovery


def test_df_each_pred() -> None:
    n_rows, n_cols = discovery.df_each_pred.shape
    assert n_rows == len(df_wbm)
    assert n_rows == len(discovery.df_preds)
    assert {*discovery.df_each_pred} == {*df_metrics}, (
        f"{discovery.df_each_pred.columns=}, expected {df_metrics.columns=}"
    )


def test_df_each_err() -> None:
    assert len(discovery.df_each_err) == len(df_wbm)
    assert len(discovery.df_each_err) == len(discovery.df_preds)
    assert {*discovery.df_each_err} == {*df_metrics, MbdKey.each_err_models}, (
        f"{discovery.df_each_err.columns=}, expected {df_metrics.columns=}"
    )
