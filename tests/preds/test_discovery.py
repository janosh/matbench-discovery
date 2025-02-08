from typing import Any

import pandas as pd
import pytest

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.preds import discovery


@pytest.fixture(autouse=True)
def mock_glob_to_df(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock glob_to_df to return test data."""
    test_data = {
        "material_id": [f"mp-{idx}" for idx in range(1, 6)],
        "e_form_pred": [0.123, 0.456, 0.789, -0.123, -0.456],
    }

    def mock_glob(*_args: Any, **_kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(test_data)

    monkeypatch.setattr("matbench_discovery.data.glob_to_df", mock_glob)


def test_df_wbm() -> None:
    for col in (MbdKey.e_form_dft, MbdKey.each_true):
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


def test_df_metrics() -> None:
    missing_cols = {*discovery.df_metrics} - {model.label for model in Model}
    assert missing_cols == set(), f"{missing_cols=}"
    assert discovery.df_metrics.T.MAE.between(0, 0.2).all(), (
        f"unexpected {discovery.df_metrics.T.MAE=}"
    )
    assert discovery.df_metrics.T.R2.between(-1.5, 1).all(), (
        f"unexpected {discovery.df_metrics.T.R2=}"
    )
    assert discovery.df_metrics.T.RMSE.between(0, 0.3).all(), (
        f"unexpected {discovery.df_metrics.T.RMSE=}"
    )
    assert discovery.df_metrics.isna().sum().sum() == 0, "NaNs in metrics"


def test_df_each_pred() -> None:
    assert len(discovery.df_each_pred) == len(df_wbm)
    assert {*discovery.df_each_pred} == {
        *discovery.df_metrics,
        MbdKey.each_mean_models,
    }, "discovery.df_each_pred has wrong columns"


def test_df_each_err() -> None:
    assert len(discovery.df_each_err) == len(df_wbm)
    assert {*discovery.df_each_err} == {
        *discovery.df_metrics,
        MbdKey.each_err_models,
    }, "discovery.df_each_err has wrong columns"
