"""Tests for package-level helpers and import-time configuration."""

import math
import os
import warnings
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Element

from matbench_discovery import ROOT, atomic_output_path, json_ready, repo_relative_path


@pytest.mark.parametrize(
    ("file_path", "expected"),
    [
        (f"{ROOT}/models/foo.json.gz", "models/foo.json.gz"),  # absolute, in repo
        ("models/foo.json.gz", "models/foo.json.gz"),  # relative, in repo
        (f"{ROOT}/a/../b.json", "b.json"),  # absolute, normalized into repo
    ],
)
def test_repo_relative_path(file_path: str, expected: str) -> None:
    """repo_relative_path returns POSIX repo-relative paths for in-repo inputs."""
    assert repo_relative_path(file_path) == expected


@pytest.mark.parametrize(
    "file_path",
    ["../outside.json", f"{ROOT}/../outside.json", "models/../../escape.json"],
    ids=["relative_escape", "absolute_outside", "relative_double_escape"],
)
def test_repo_relative_path_rejects_outside_root(file_path: str) -> None:
    """Paths escaping the repo root are rejected, including relative ones."""
    with pytest.raises(ValueError, match="must be inside repo root"):
        repo_relative_path(file_path)


def test_atomic_output_path_keeps_old_file_on_error(tmp_path: Path) -> None:
    """A failed write keeps the previous file and leaves no temp file behind."""
    target = f"{tmp_path}/out.txt"
    with open(target, mode="w", encoding="utf-8") as file:
        file.write("old")

    def write_then_fail() -> None:
        """Write the new content, then crash before the rename."""
        with atomic_output_path(target) as temporary_path:
            with open(temporary_path, mode="w", encoding="utf-8") as file:
                file.write("new")
            raise RuntimeError

    with pytest.raises(RuntimeError):
        write_then_fail()
    with open(target, encoding="utf-8") as file:
        assert file.read() == "old"
    assert os.listdir(tmp_path) == ["out.txt"]


@pytest.mark.parametrize("nan_to_none", [True, False])
def test_json_ready(nan_to_none: bool) -> None:
    """NumPy, tuple and key conversion; non-finite floats become None or raise."""
    data = {1: (np.int64(2), np.array([0.5, 1.5])), "flag": np.bool_(1)}
    assert json_ready(data, nan_to_none=nan_to_none) == {
        "1": [2, [0.5, 1.5]],
        "flag": True,
    }
    if nan_to_none:
        assert json_ready([np.float32("nan"), math.inf]) == [None, None]
    else:
        with pytest.raises(ValueError, match="not JSON compliant: nan"):
            json_ready({"x": [np.nan]}, nan_to_none=False)
    with pytest.raises(TypeError, match="object is not JSON serializable"):
        json_ready(object(), nan_to_none=nan_to_none)


def test_pymatgen_electronegativity_warning_filtered() -> None:
    """Import-time filters silence pymatgen's missing-electronegativity warning."""
    # call the uncached getter: Element.X caches after its first (warning) access
    with warnings.catch_warnings(record=True) as caught:
        type(Element.Ne).X.func(Element.Ne)
    assert not [warn for warn in caught if "electronegativity" in str(warn.message)]
