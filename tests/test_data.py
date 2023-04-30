from __future__ import annotations

import json
import os
from pathlib import Path
from random import random
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pytest import CaptureFixture

from matbench_discovery import FIGSHARE, ROOT
from matbench_discovery.data import (
    DATA_FILES,
    as_dict_handler,
    df_wbm,
    figshare_versions,
    glob_to_df,
    load_train_test,
)

with open(f"{FIGSHARE}/{figshare_versions[-1]}.json") as file:
    figshare_urls = json.load(file)

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


@pytest.mark.parametrize(
    "data_names, hydrate",
    [
        (["wbm_summary"], True),
        (["wbm_initial_structures"], True),
        (["wbm_computed_structure_entries"], False),
        (["wbm_summary", "wbm_initial_structures"], True),
        (["mp_elemental_ref_entries"], True),
        (["mp_energies"], True),
    ],
)
def test_load_train_test(
    data_names: list[str],
    hydrate: bool,
    dummy_df_serialized: pd.DataFrame,
    capsys: CaptureFixture[str],
    tmp_path: Path,
) -> None:
    # intercept HTTP requests to GitHub raw user content and return dummy df instead
    with patch("matbench_discovery.data.pd.read_csv") as read_csv, patch(
        "matbench_discovery.data.pd.read_json"
    ) as read_json:
        # dummy df with Structures and ComputedStructureEntries
        read_json.return_value = dummy_df_serialized
        # dummy df with random floats and material_id column
        read_csv.return_value = pd._testing.makeDataFrame().reset_index(
            names="material_id"
        )
        out = load_train_test(
            data_names,
            hydrate=hydrate,
            # test both str and Path for cache_dir
            cache_dir=str(tmp_path) if random() < 0.5 else tmp_path,
        )

    stdout, _stderr = capsys.readouterr()

    expected_outs = [
        f"Downloading {key!r} from {figshare_urls[key]}" for key in data_names
    ]
    for expected_out in expected_outs:
        assert expected_out in stdout

    # check we called read_csv/read_json once for each data_name
    assert read_json.call_count + read_csv.call_count == len(data_names)

    if len(data_names) > 1:
        assert isinstance(out, dict)
        assert list(out) == data_names
        for key, df in out.items():
            assert isinstance(df, pd.DataFrame), f"{key} not a DataFrame but {type(df)}"
    else:
        assert isinstance(out, pd.DataFrame), f"{data_names[0]} not a DataFrame"

    # test that df loaded from cache is the same as initial df
    from_cache = load_train_test(data_names, hydrate=hydrate, cache_dir=tmp_path)
    if len(data_names) > 1:
        for key, df in from_cache.items():
            pd.testing.assert_frame_equal(df, out[key])
    else:
        pd.testing.assert_frame_equal(out, from_cache)


def test_load_train_test_raises(tmp_path: Path) -> None:
    # bad data name
    with pytest.raises(ValueError, match=f"must be subset of {set(DATA_FILES)}"):
        load_train_test(["bad-data-name"])

    # bad_version
    version = "invalid-version"
    with pytest.raises(ValueError) as exc_info:
        load_train_test("wbm_summary", version=version, cache_dir=tmp_path)

    assert (
        str(exc_info.value)
        == f"Unexpected {version=}. Must be one of {figshare_versions}."
    )
    assert os.listdir(tmp_path) == [], "cache_dir should be empty"


def test_load_train_test_doc_str() -> None:
    doc_str = load_train_test.__doc__
    assert isinstance(doc_str, str)  # mypy type narrowing

    # check that we link to the right data description page
    with open(f"{ROOT}/site/package.json") as file:
        pkg = json.load(file)  # get repo URL from package.json
    route = "/contribute"
    assert f"{pkg['homepage']}/contribute" in doc_str
    assert os.path.isdir(f"{ROOT}/site/src/routes/{route}")


wbm_summary_expected_cols = {
    "bandgap_pbe",
    "e_form_per_atom_mp2020_corrected",
    "e_form_per_atom_uncorrected",
    "e_form_per_atom_wbm",
    "e_above_hull_wbm",
    "formula",
    "n_sites",
    "uncorrected_energy",
    "uncorrected_energy_from_cse",
    "volume",
    "wyckoff_spglib",
}


# TODO skip this test if offline
# @pytest.mark.skipif(online, reason="requires internet connection")
@pytest.mark.parametrize(
    "file_key, version, expected_shape, expected_cols",
    [
        ("mp_elemental_ref_entries", figshare_versions[-1], (5, 89), set()),
        pytest.param(
            "wbm_summary",
            figshare_versions[-1],
            (256963, 15),
            wbm_summary_expected_cols,
            marks=pytest.mark.slow,  # run pytest -m 'slow' to select this marker
        ),
        pytest.param(
            # large file but needed to test loading compressed JSON from URL
            "mp_computed_structure_entries",
            figshare_versions[-1],
            (154718, 1),
            {"entry"},
            marks=pytest.mark.slow,
        ),
    ],
)
def test_load_train_test_no_mock(
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
    df = load_train_test(file_key, version=version, cache_dir=tmp_path)
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
        f"Downloading {file_key!r} from {figshare_urls[file_key]}\nCached "
        f"{file_key!r} to {cache_path!r}" in stdout
    )

    # test that df loaded from cache is the same as initial df
    pd.testing.assert_frame_equal(
        df, load_train_test(file_key, version=version, cache_dir=tmp_path)
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
    assert df_wbm.shape == (256963, 16)
    assert df_wbm.index.name == "material_id"
    assert set(df_wbm) > {"bandgap_pbe", "formula", "material_id"}


@pytest.mark.parametrize("pattern", ["tmp/*df.csv", "tmp/*df.json"])
def test_glob_to_df(pattern: str) -> None:
    try:
        df = pd._testing.makeMixedDataFrame()

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
