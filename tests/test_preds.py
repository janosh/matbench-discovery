from matbench_discovery.data import PRED_FILES
from matbench_discovery.preds import (
    df_each_err,
    df_each_pred,
    df_metrics,
    df_wbm,
    e_form_col,
    each_true_col,
)


def test_df_wbm() -> None:
    for col in [e_form_col, each_true_col]:
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


def test_df_metrics() -> None:
    assert {*df_metrics} == {*PRED_FILES}
    assert df_metrics.T.MAE.between(0, 0.2).all(), f"unexpected {df_metrics.T.MAE=}"
    assert df_metrics.T.R2.between(-0.65, 1).all(), f"unexpected {df_metrics.T.R2=}"
    assert df_metrics.T.RMSE.between(0, 0.25).all(), f"unexpected {df_metrics.T.RMSE=}"
    assert df_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_df_each_pred() -> None:
    assert len(df_each_pred) == len(df_wbm)
    assert {*df_each_pred} == {*df_metrics}, "df_each_pred has wrong columns"
    assert all(df_each_pred.isna().mean() < 0.05), "too many NaNs in df_each_pred"


def test_df_each_err() -> None:
    assert len(df_each_err) == len(df_wbm)
    assert {*df_each_err} == {*df_metrics}, "df_each_err has wrong columns"
    assert all(df_each_err.isna().mean() < 0.05), "too many NaNs in df_each_err"
