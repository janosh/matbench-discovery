from __future__ import annotations

import pytest

from matbench_discovery.plot_scripts import (
    data_paths,
    df_wbm,
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
    assert len(data_paths) >= 9
    assert all(path.startswith(("models/", "data/")) for path in data_paths.values())
