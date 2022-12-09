from __future__ import annotations

import os

import pandas as pd
import pytest

from matbench_discovery import ROOT
from matbench_discovery.load_preds import (
    DATA_PATHS,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
)


def test_df_wbm() -> None:
    assert len(df_wbm) == 256963
    assert df_wbm.shape >= (256963, 18)
    assert df_wbm.index.name == "material_id"
    assert set(df_wbm) > {"bandgap_pbe", "formula", "material_id"}


@pytest.mark.parametrize("models", [[], ["Wrenformer"]])
def test_load_df_wbm_with_preds(models: list[str]) -> None:
    df = load_df_wbm_with_preds(models=models)
    assert len(df) == len(df_wbm)
    assert list(df) == list(df_wbm) + models + [f"{model}_std" for model in models]
    assert df.index.name == "material_id"

    for model_name in models:
        assert model_name in df
        assert df[model_name].isna().sum() == 0


def test_load_df_wbm_with_preds_raises() -> None:
    with pytest.raises(ValueError, match="Unknown models: foo"):
        load_df_wbm_with_preds(models=["foo"])


def test_data_paths() -> None:
    assert len(DATA_PATHS) >= 6
    assert all(path.startswith(("models/", "data/")) for path in DATA_PATHS.values())


@pytest.mark.parametrize("pattern", ["tmp/*df.csv", "tmp/*df.json"])
def test_glob_to_df(pattern: str) -> None:
    try:
        df = pd.util.testing.makeMixedDataFrame()

        os.makedirs(f"{ROOT}/tmp", exist_ok=True)
        df.to_csv(f"{ROOT}/tmp/dummy_df.csv", index=False)
        df.to_json(f"{ROOT}/tmp/dummy_df.json")

        df_out = glob_to_df(pattern)
        assert df_out.shape == df.shape
        assert list(df_out) == list(df)

        with pytest.raises(FileNotFoundError):
            glob_to_df("foo")
    finally:
        os.remove(f"{ROOT}/tmp/dummy_df.csv")
        os.remove(f"{ROOT}/tmp/dummy_df.json")
