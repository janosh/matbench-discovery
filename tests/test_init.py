"""Tests for the repo_relative_path package helper."""

import pytest

from matbench_discovery import ROOT, repo_relative_path


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
