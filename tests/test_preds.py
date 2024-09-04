import os

import numpy as np
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.preds import (
    Model,
    df_each_err,
    df_each_pred,
    df_metrics,
    load_df_wbm_with_preds,
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


@pytest.mark.parametrize("models", [[], ["Wrenformer"]])
@pytest.mark.parametrize("max_error_threshold", [None, 5.0, 1.0])
def test_load_df_wbm_with_preds(
    models: list[str], max_error_threshold: float | None
) -> None:
    df_wbm_with_preds = load_df_wbm_with_preds(
        models=models, max_error_threshold=max_error_threshold
    )
    assert len(df_wbm_with_preds) == len(df_wbm)

    assert list(df_wbm_with_preds) == list(df_wbm) + [
        Model.label_map.get(model, model) for model in models
    ]
    assert df_wbm_with_preds.index.name == Key.mat_id

    for model_name in models:
        assert model_name in df_wbm_with_preds
        if max_error_threshold is not None:
            # Check if predictions exceeding the threshold are filtered out
            error = abs(
                df_wbm_with_preds[model_name] - df_wbm_with_preds[MbdKey.e_form_dft]
            )
            assert np.all(error[~error.isna()] <= max_error_threshold)
        else:
            # If no threshold is set, all predictions should be present
            assert df_wbm_with_preds[model_name].isna().sum() == 0


def test_load_df_wbm_max_error_threshold() -> None:
    models = {Model.mace.label: 38}  # num missing preds for default max_error_threshold
    df_no_thresh = load_df_wbm_with_preds(models=list(models))
    df_high_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=10)
    df_low_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=0.1)

    for model, n_missing in models.items():
        assert df_no_thresh[model].isna().sum() == n_missing
        assert df_high_thresh[model].isna().sum() <= df_no_thresh[model].isna().sum()
        assert df_high_thresh[model].isna().sum() <= df_low_thresh[model].isna().sum()


def test_load_df_wbm_with_preds_raises() -> None:
    with pytest.raises(ValueError, match="unknown_models='foo'"):
        load_df_wbm_with_preds(models=["foo"])

    with pytest.raises(
        ValueError, match="max_error_threshold must be a positive number"
    ):
        load_df_wbm_with_preds(max_error_threshold=-1.0)


def test_pred_files() -> None:
    assert len(Model) >= 6
    assert all(
        file.path.endswith((".csv", ".csv.gz", ".json", ".json.gz")) for file in Model
    )
    for model in Model:
        path = model.path
        msg = f"Missing preds file for {model=}, expected at {path=}"
        assert os.path.isfile(path), msg
