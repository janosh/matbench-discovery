from __future__ import annotations

import os
from tempfile import TemporaryDirectory
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure

from matbench_discovery import ROOT
from matbench_discovery.data import (
    DATA_FILENAMES,
    PRED_FILENAMES,
    RAW_REPO_URL,
    as_dict_handler,
    chunks,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
    load_train_test,
)

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


@pytest.mark.parametrize(
    "parts, cache_dir, hydrate",
    [
        (["wbm-summary"], None, True),
        (["wbm-initial-structures"], TemporaryDirectory().name, True),
        (["wbm-computed-structure-entries"], None, False),
        (["wbm-summary", "wbm-initial-structures"], TemporaryDirectory().name, True),
        (["mp-elemental-ref-energies"], None, True),
        (["mp-energies"], None, True),
    ],
)
def test_load_train_test(
    parts: list[str],
    cache_dir: str | None,
    hydrate: bool,
    dummy_df_with_structures: pd.DataFrame,
    capsys: pytest.CaptureFixture,
) -> None:
    # intercept HTTP requests to GitHub raw user content and return dummy df instead
    with patch("matbench_discovery.data.pd.read_csv") as read_csv, patch(
        "matbench_discovery.data.pd.read_json"
    ) as read_json:
        read_csv.return_value = read_json.return_value = dummy_df_with_structures
        out = load_train_test(parts, cache_dir=cache_dir, hydrate=hydrate)

    stdout, stderr = capsys.readouterr()

    assert (
        "\n".join(
            f"Downloading {part} from {RAW_REPO_URL}/1.0.0/data/{DATA_FILENAMES[part]}"
            for part in parts
        )
        in stdout
    )
    assert "" == stderr

    assert read_json.call_count + read_csv.call_count == len(parts)

    if len(parts) > 1:
        assert isinstance(out, dict)
        assert list(out) == parts
        for df in out.values():
            assert isinstance(df, pd.DataFrame)
    else:
        assert isinstance(out, pd.DataFrame)


def test_load_train_test_raises() -> None:
    with pytest.raises(
        ValueError,
        match=f"must be subset of {set(DATA_FILENAMES)}",
    ):
        load_train_test(["invalid-part"])

    with pytest.raises(
        ValueError, match="Only version 1 currently available, got version=2"
    ):
        load_train_test(version=2)


def test_load_train_test_doc_str() -> None:
    doc_str = load_train_test.__doc__
    assert isinstance(doc_str, str)  # mypy type narrowing

    assert all(key in doc_str for key in DATA_FILENAMES)

    # TODO refactor to load site URL from site/package.json for SSoT
    assert "https://matbench-discovery.janosh.dev" in doc_str


def test_chunks() -> None:
    assert list(chunks([], 1)) == []
    assert list(chunks([1], 1)) == [[1]]
    assert list(chunks([1, 2], 1)) == [[1], [2]]
    assert list(chunks([1, 2, 3], 1)) == [[1], [2], [3]]
    assert list(chunks([1, 2, 3], 2)) == [[1, 2], [3]]
    assert list(chunks(range(1, 4), 2)) == [range(1, 3), range(3, 4)]
    assert list(chunks(range(1, 5), 2)) == [range(1, 3), range(3, 5)]
    assert list(chunks(range(1, 5), 3)) == [range(1, 4), range(4, 5)]


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


def test_pred_filenames() -> None:
    assert len(PRED_FILENAMES) >= 6
    assert all(
        path.startswith(("models/", "data/")) for path in PRED_FILENAMES.values()
    )


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
