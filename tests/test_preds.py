from matbench_discovery.data import PRED_FILENAMES
from matbench_discovery.preds import (
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
    assert {*df_metrics} == {*PRED_FILENAMES}
    assert df_metrics.T.MAE.between(0, 0.2).all(), f"unexpected {df_metrics.T.MAE=}"
    assert df_metrics.T.R2.between(-0.65, 1).all(), f"unexpected {df_metrics.T.R2=}"
    assert df_metrics.T.RMSE.between(0, 0.25).all(), f"unexpected {df_metrics.T.RMSE=}"
    assert df_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_df_each_pred() -> None:
    assert len(df_each_pred) == len(df_wbm)
    assert (
        {*df_each_pred} == {*df_metrics} < {*df_wbm}
    ), "df_each_pred has wrong columns"
    assert all(
        df_each_pred.isna().sum() / len(df_each_pred) < 0.05
    ), "too many NaNs in df_each_pred"
