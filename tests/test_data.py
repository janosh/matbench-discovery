from __future__ import annotations

import json
import os
from random import random
from typing import TYPE_CHECKING, Any
from unittest.mock import patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pytest import CaptureFixture

from matbench_discovery import FIGSHARE_DIR, ROOT, Key
from matbench_discovery.data import (
    DATA_FILES,
    as_dict_handler,
    df_wbm,
    figshare_versions,
    glob_to_df,
    load,
)

if TYPE_CHECKING:
    from pathlib import Path


with open(f"{FIGSHARE_DIR}/{figshare_versions[-1]}.json") as file:
    figshare_urls = json.load(file)["files"]

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


@pytest.mark.parametrize(
    "key, hydrate",
    [
        ("wbm_summary", True),
        ("wbm_initial_structures", True),
        ("wbm_computed_structure_entries", False),
        ("mp_elemental_ref_entries", True),
        ("mp_energies", True),
    ],
)
def test_load(
    df_float: pd.DataFrame,
    # df with Structures and ComputedStructureEntries as dicts
    df_with_pmg_objects: pd.DataFrame,
    capsys: CaptureFixture[str],
    tmp_path: Path,
    key: str,
    hydrate: bool,
) -> None:
    filepath = DATA_FILES[key]
    # intercept HTTP requests and write dummy df to disk instead
    with patch("urllib.request.urlretrieve") as url_retrieve:
        writer = df_with_pmg_objects.to_json if ".json" in filepath else df_float.to_csv
        url_retrieve.side_effect = lambda _url, path: writer(path)
        out = load(
            key,
            hydrate=hydrate,
            # test both str and Path for cache_dir
            cache_dir=str(tmp_path) if random() < 0.5 else tmp_path,
        )

    stdout, _stderr = capsys.readouterr()

    assert f"Downloading {key!r} from {figshare_urls[key][0]}" in stdout

    # check we called read_csv/read_json once for each data_name
    assert url_retrieve.call_count == 1

    assert isinstance(out, pd.DataFrame), f"{key} not a DataFrame"

    # test that df loaded from cache is the same as initial df
    from_cache = load(key, hydrate=hydrate, cache_dir=tmp_path)
    pd.testing.assert_frame_equal(out, from_cache)


def test_load_raises(tmp_path: Path) -> None:
    key = "bad-key"
    with pytest.raises(ValueError) as exc:  # noqa: PT011
        load(key)

    assert f"Unknown {key=}, must be one of {list(DATA_FILES)}" in str(exc.value)

    version = "invalid-version"
    with pytest.raises(ValueError) as exc:  # noqa: PT011
        load("wbm_summary", version=version, cache_dir=tmp_path)

    assert (
        str(exc.value) == f"Unexpected {version=}. Must be one of {figshare_versions}."
    )
    assert os.listdir(tmp_path) == [], "cache_dir should be empty"


def test_load_doc_str() -> None:
    doc_str = load.__doc__
    assert isinstance(doc_str, str)  # mypy type narrowing

    # check that we link to the right data description page
    with open(f"{ROOT}/site/package.json") as file:
        pkg = json.load(file)  # get repo URL from package.json
    route = "/contribute"
    assert f"{pkg['homepage']}/contribute" in doc_str
    assert os.path.isdir(f"{ROOT}/site/src/routes/{route}")


wbm_summary_expected_cols = {
    "uncorrected_energy_from_cse",
    Key.bandgap_pbe,
    Key.dft_energy,
    Key.e_form_raw,
    Key.e_form_wbm,
    Key.e_form,
    Key.each_wbm,
    Key.formula,
    Key.n_sites,
    Key.volume,
    Key.wyckoff,
}


# TODO skip this test if offline
# @pytest.mark.skipif(online, reason="requires internet connection")
@pytest.mark.parametrize(
    "file_key, version, expected_shape, expected_cols",
    [
        ("mp_elemental_ref_entries", figshare_versions[-1], (9, 89), set()),
        pytest.param(
            "wbm_summary",
            figshare_versions[-1],
            (256_963, 17),
            wbm_summary_expected_cols,
            marks=pytest.mark.slow,  # run pytest -m 'slow' to select this marker
        ),
        pytest.param(
            # large file but needed to test loading compressed JSON from URL
            "mp_computed_structure_entries",
            figshare_versions[-1],
            (154_718, 1),
            {"entry"},
            marks=pytest.mark.very_slow,
        ),
    ],
)
def test_load_no_mock(
    file_key: str,
    version: str,
    expected_shape: tuple[int, int],
    expected_cols: set[str],
    capsys: CaptureFixture[str],
    tmp_path: Path,
) -> None:
    assert os.listdir(tmp_path) == [], "cache_dir should be empty"
    # This function runs the download from Figshare for real hence takes some time and
    # requires being online
    df = load(file_key, version=version, cache_dir=tmp_path)
    assert len(os.listdir(tmp_path)) == 1, "cache_dir should have one file"
    assert df.shape == expected_shape
    assert (
        set(df) >= expected_cols
    ), f"Loaded df missing columns {expected_cols - set(df)}"

    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    rel_path = getattr(type(DATA_FILES), file_key)
    cache_path = f"{tmp_path}/{rel_path}"
    assert (
        f"Downloading {file_key!r} from {figshare_urls[file_key][0]}\nCached "
        f"{file_key!r} to {cache_path!r}" in stdout
    )

    # test that df loaded from cache is the same as initial df
    pd.testing.assert_frame_equal(
        df, load(file_key, version=version, cache_dir=tmp_path)
    )

    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    assert stdout == f"Loading {file_key!r} from cached file at {cache_path!r}\n"


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
