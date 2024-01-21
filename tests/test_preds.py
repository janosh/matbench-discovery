import os

import pytest

from matbench_discovery import Key
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import (
    PRED_FILES,
    df_each_err,
    df_each_pred,
    df_metrics,
    load_df_wbm_with_preds,
)


def test_df_wbm() -> None:
    for col in (Key.e_form, Key.each_true):
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


def test_df_metrics() -> None:
    assert {*df_metrics} >= {*PRED_FILES}
    assert df_metrics.T.MAE.between(0, 0.2).all(), f"unexpected {df_metrics.T.MAE=}"
    assert df_metrics.T.R2.between(-1.5, 1).all(), f"unexpected {df_metrics.T.R2=}"
    assert df_metrics.T.RMSE.between(0, 0.3).all(), f"unexpected {df_metrics.T.RMSE=}"
    assert df_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_df_each_pred() -> None:
    assert len(df_each_pred) == len(df_wbm)
    assert {*df_each_pred} == {
        *df_metrics,
        Key.model_mean_each,
    }, "df_each_pred has wrong columns"
    assert all(df_each_pred.isna().mean() < 0.05), "too many NaNs in df_each_pred"


def test_df_each_err() -> None:
    assert len(df_each_err) == len(df_wbm)
    assert {*df_each_err} == {
        *df_metrics,
        Key.model_mean_err,
    }, "df_each_err has wrong columns"
    assert all(df_each_err.isna().mean() < 0.05), "too many NaNs in df_each_err"


@pytest.mark.parametrize("models", [[], ["Wrenformer"]])
def test_load_df_wbm_with_preds(models: list[str]) -> None:
    df = load_df_wbm_with_preds(models)
    assert len(df) == len(df_wbm)
    assert list(df) == list(df_wbm) + models + [f"{model}_std" for model in models]
    assert df.index.name == Key.mat_id

    for model_name in models:
        assert model_name in df
        assert df[model_name].isna().sum() == 0


def test_load_df_wbm_with_preds_raises() -> None:
    with pytest.raises(ValueError, match="Unknown models: foo"):
        load_df_wbm_with_preds(models=["foo"])


def test_pred_files() -> None:
    assert len(PRED_FILES) >= 6
    assert all(
        path.endswith((".csv", ".csv.gz", ".json", ".json.gz"))
        for path in PRED_FILES.values()
    )
    for model, path in PRED_FILES.items():
        msg = f"Missing preds file for {model=}, expected at {path=}"
        assert os.path.isfile(path), msg
