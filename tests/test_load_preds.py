from __future__ import annotations

import os

import pandas as pd
import pytest

from matbench_discovery import ROOT
from matbench_discovery.load_preds import (
    data_paths,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
    load_model_preds,
)


def test_df_wbm() -> None:
    assert len(df_wbm) == 257487
    assert df_wbm.shape >= (257487, 18)
    assert df_wbm.index.name == "material_id"
    assert set(df_wbm) > {"bandgap_pbe", "formula", "material_id"}


def test_load_df_wbm_with_preds() -> None:
    df = load_df_wbm_with_preds(models=[])
    assert df.shape == df_wbm.shape
    assert df.index.name == "material_id"


def test_load_model_preds() -> None:
    dfs = load_model_preds(models=[])
    assert dfs == {}

    with pytest.raises(ValueError):
        load_model_preds(models=["foo"])


def test_data_paths() -> None:
    assert len(data_paths) >= 8
    assert all(path.startswith(("models/", "data/")) for path in data_paths.values())


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