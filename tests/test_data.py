import json
import os
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatviz.enums import Key

from matbench_discovery import FIGSHARE_DIR
from matbench_discovery.data import (
    as_dict_handler,
    df_wbm,
    figshare_versions,
    glob_to_df,
)

with open(f"{FIGSHARE_DIR}/{figshare_versions[-1]}.json") as file:
    figshare_urls = json.load(file)["files"]

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


def test_as_dict_handler() -> None:
    class C:
        def as_dict(self) -> dict[str, Any]:
            return {"foo": "bar"}

    assert as_dict_handler(C()) == {"foo": "bar"}
    assert as_dict_handler(1) is None
    assert as_dict_handler("foo") is None
    assert as_dict_handler([1, 2, 3]) is None
    assert as_dict_handler({"foo": "bar"}) is None


def test_df_wbm() -> None:
    assert df_wbm.shape == (256_963, 18)
    assert df_wbm.index.name == Key.mat_id
    assert set(df_wbm) > {Key.formula, Key.mat_id, Key.bandgap_pbe}


@pytest.mark.parametrize("pattern", ["*df.csv", "*df.json"])
def test_glob_to_df(pattern: str, tmp_path: Path, df_mixed: pd.DataFrame) -> None:
    os.makedirs(f"{tmp_path}", exist_ok=True)
    df_mixed.to_csv(f"{tmp_path}/dummy_df.csv", index=False)
    df_mixed.to_json(f"{tmp_path}/dummy_df.json")

    df_out = glob_to_df(f"{tmp_path}/{pattern}")
    assert df_out.shape == df_mixed.shape
    assert list(df_out) == list(df_mixed)

    with pytest.raises(FileNotFoundError):
        glob_to_df("foo")
