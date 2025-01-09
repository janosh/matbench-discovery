import os

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.preds.discovery import (
    Model,
    df_each_err,
    df_each_pred,
    df_metrics,
)


def test_df_wbm() -> None:
    for col in (MbdKey.e_form_dft, MbdKey.each_true):
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


def test_df_metrics() -> None:
    missing_cols = {*df_metrics} - {model.label for model in Model}
    assert missing_cols == set(), f"{missing_cols=}"
    assert df_metrics.T.MAE.between(0, 0.2).all(), f"unexpected {df_metrics.T.MAE=}"
    assert df_metrics.T.R2.between(-1.5, 1).all(), f"unexpected {df_metrics.T.R2=}"
    assert df_metrics.T.RMSE.between(0, 0.3).all(), f"unexpected {df_metrics.T.RMSE=}"
    assert df_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_df_each_pred() -> None:
    assert len(df_each_pred) == len(df_wbm)
    assert {*df_each_pred} == {
        *df_metrics,
        MbdKey.each_mean_models,
    }, "df_each_pred has wrong columns"
    assert all(df_each_pred.isna().mean() < 0.05), "too many NaNs in df_each_pred"


def test_df_each_err() -> None:
    assert len(df_each_err) == len(df_wbm)
    assert {*df_each_err} == {
        *df_metrics,
        MbdKey.each_err_models,
    }, "df_each_err has wrong columns"
    # get model names that have more than 5 percent NaNs
    models_with_too_many_nans = df_each_err.columns[
        df_each_err.isna().mean() > 0.05
    ].tolist()
    assert len(models_with_too_many_nans) == 0, (
        "Some models have too many NaNs in df_each_err:\n"
        + "\n".join(
            f"- {model}: {df_each_err[model].isna().mean():.2%}"
            for model in models_with_too_many_nans
        )
    )


def test_pred_files() -> None:
    assert len(Model) >= 6
    for model in Model:
        pred_path = model.discovery_path
        assert pred_path.endswith(".csv.gz")
        assert os.path.isfile(pred_path), (
            f"discovery pred file for {model=} not found, expected at {pred_path}"
        )
