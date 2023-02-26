from __future__ import annotations

import os
import urllib.request
from pathlib import Path
from random import random
from tempfile import TemporaryDirectory
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pytest import CaptureFixture

from matbench_discovery import ROOT
from matbench_discovery.data import (
    DATA_FILENAMES,
    PRED_FILENAMES,
    RAW_REPO_URL,
    as_dict_handler,
    df_wbm,
    glob_to_df,
    load_df_wbm_preds,
    load_train_test,
)

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)

try:
    website_down = urllib.request.urlopen(RAW_REPO_URL).status != 200
except Exception:
    website_down = True


@pytest.mark.parametrize(
    "data_names, hydrate",
    [
        (["wbm-summary"], True),
        (["wbm-initial-structures"], True),
        (["wbm-computed-structure-entries"], False),
        (["wbm-summary", "wbm-initial-structures"], True),
        (["mp-elemental-ref-energies"], True),
        (["mp-energies"], True),
    ],
)
def test_load_train_test(
    data_names: list[str],
    hydrate: bool,
    dummy_df_with_structures: pd.DataFrame,
    capsys: CaptureFixture[str],
    tmp_path: Path,
) -> None:
    # intercept HTTP requests to GitHub raw user content and return dummy df instead
    with patch("matbench_discovery.data.pd.read_csv") as read_csv, patch(
        "matbench_discovery.data.pd.read_json"
    ) as read_json:
        read_csv.return_value = read_json.return_value = dummy_df_with_structures
        out = load_train_test(
            data_names,
            hydrate=hydrate,
            # test both str and Path cache_dir
            cache_dir=TemporaryDirectory().name if random() < 0.5 else tmp_path,
        )

    stdout, stderr = capsys.readouterr()

    expected_out = "\n".join(
        f"Downloading {name!r} from {RAW_REPO_URL}/1.0.0/data/{DATA_FILENAMES[name]}"
        for name in data_names
    )
    assert expected_out in stdout
    assert stderr == ""

    assert read_json.call_count + read_csv.call_count == len(data_names)

    if len(data_names) > 1:
        assert isinstance(out, dict)
        assert list(out) == data_names
        for df in out.values():
            assert isinstance(df, pd.DataFrame)
    else:
        assert isinstance(out, pd.DataFrame)


def test_load_train_test_raises(tmp_path: Path) -> None:
    # bad data name
    with pytest.raises(ValueError, match=f"must be subset of {set(DATA_FILENAMES)}"):
        load_train_test(["bad-data-name"])

    # bad_version
    version = "not-a-real-branch"
    with pytest.raises(ValueError) as exc_info:
        load_train_test("wbm-summary", version=version, cache_dir=tmp_path)

    assert (
        str(exc_info.value)
        == "Bad url='https://raw.githubusercontent.com/janosh/matbench-discovery"
        f"/{version}/data/wbm/2022-10-19-wbm-summary.csv'"
    )


def test_load_train_test_doc_str() -> None:
    doc_str = load_train_test.__doc__
    assert isinstance(doc_str, str)  # mypy type narrowing

    for name in DATA_FILENAMES:
        assert name in doc_str, f"Missing data {name=} in load_train_test() docstring"

    route = "/contribute"
    assert os.path.isdir(f"{ROOT}/site/src/routes/{route}")


@pytest.mark.skipif(website_down, reason=f"{RAW_REPO_URL} unreachable")
@pytest.mark.parametrize("version", ["main"])  # , "d00d475"
def test_load_train_test_no_mock(
    version: str, capsys: CaptureFixture[str], tmp_path: Path
) -> None:
    # this function runs the download from GitHub raw user content for real
    # hence takes some time and requires being online
    df_wbm = load_train_test("wbm-summary", version=version, cache_dir=tmp_path)
    assert df_wbm.shape == (256963, 14)
    assert set(df_wbm) > {
        "bandgap_pbe",
        "e_form_per_atom_mp2020_corrected",
        "e_form_per_atom_uncorrected",
        "e_form_per_atom_wbm",
        "e_hull_wbm",
        "formula",
        "n_sites",
        "uncorrected_energy",
        "uncorrected_energy_from_cse",
        "volume",
        "wyckoff_spglib",
    }, "Loaded df missing columns"

    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    assert (
        stdout
        == "Downloading 'wbm-summary' from https://raw.githubusercontent.com/janosh"
        f"/matbench-discovery/{version}/data/wbm/2022-10-19-wbm-summary.csv\n"
    )

    df_wbm = load_train_test("wbm-summary", version=version, cache_dir=tmp_path)

    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    assert (
        stdout
        == f"Loading 'wbm-summary' from cached file at '{tmp_path}/main/wbm/2022-10-19-"
        "wbm-summary.csv'\n"
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
    assert df_wbm.shape == (256963, 14)
    assert df_wbm.index.name == "material_id"
    assert set(df_wbm) > {"bandgap_pbe", "formula", "material_id"}


@pytest.mark.parametrize("models", [[], ["Wrenformer"]])
def test_load_df_wbm_with_preds(models: list[str]) -> None:
    df = load_df_wbm_preds(models)
    assert len(df) == len(df_wbm)
    assert list(df) == list(df_wbm) + models + [f"{model}_std" for model in models]
    assert df.index.name == "material_id"

    for model_name in models:
        assert model_name in df
        assert df[model_name].isna().sum() == 0


def test_load_df_wbm_with_preds_raises() -> None:
    with pytest.raises(ValueError, match="Unknown models: foo"):
        load_df_wbm_preds(models=["foo"])


def test_pred_filenames() -> None:
    assert len(PRED_FILENAMES) >= 6
    assert all(path.endswith((".csv", ".json")) for path in PRED_FILENAMES.values())


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
