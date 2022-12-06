from __future__ import annotations

from tempfile import TemporaryDirectory
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure

from matbench_discovery.data import (
    DATA_FILENAMES,
    RAW_REPO_URL,
    as_dict_handler,
    chunks,
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
def test_load_wbm(
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


def test_load_wbm_raises() -> None:
    with pytest.raises(
        ValueError,
        match=f"must be subset of {set(DATA_FILENAMES)}",
    ):
        load_train_test(["invalid-part"])

    with pytest.raises(
        ValueError, match="Only version 1 currently available, got version=2"
    ):
        load_train_test(version=2)


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
