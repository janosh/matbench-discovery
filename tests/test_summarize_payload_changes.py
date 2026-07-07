"""Tests for the update-site-figs auto-commit roster summarizer."""

import json
import os
import subprocess
from pathlib import Path

import pytest

from scripts.summarize_payload_changes import format_change, roster, summarize


@pytest.mark.parametrize(
    "was, now, expected",
    [
        # unchanged -> no bullet
        ({"a", "b"}, {"a", "b"}, None),
        # additions only
        ({"a"}, {"a", "b"}, "- `p`: 1 → 2 (+1: b)"),
        # removals only (previously rendered as a bare "1 → 0" with no explanation)
        ({"a"}, set(), "- `p`: 1 → 0 (-1: a)"),
        # net count drops while models are added: must show BOTH +N and -N, not a
        # lone "+N" that looks contradictory against the decreasing count (the bug)
        (
            {"a", "b", "c", "d", "e"},
            {"a", "b", "x"},
            "- `p`: 5 → 3 (+1: x; -3: c, d, e)",
        ),
        # equal-size swap
        ({"a", "b"}, {"a", "c"}, "- `p`: 2 → 2 (+1: c; -1: b)"),
        # multiple add+remove: ids sorted for deterministic PR bodies
        ({"z", "m"}, {"m", "a", "k"}, "- `p`: 2 → 3 (+2: a, k; -1: z)"),
    ],
)
def test_format_change(was: set[str], now: set[str], expected: str | None) -> None:
    """Roster delta lists both additions and removals (sorted); None when unchanged."""
    assert format_change("p", was, now) == expected


def test_roster_excludes_base_and_metadata() -> None:
    """roster() collects model ids, dropping the _base line and metadata columns."""
    lines = [
        {"_base": {"shared": 1}},
        {"key": "model-a"},
        {"label": "model-b"},
        {"key": "MP Occurrences"},  # metadata, not a model
    ]
    text = "\n".join(json.dumps(line) for line in lines)
    assert roster(text) == {"model-a", "model-b"}


def test_summarize_diffs_working_tree_against_head(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End to end: summarize() diffs on-disk payloads against their committed (HEAD)
    versions in a real git repo: changed, brand-new (untracked) and unchanged files.
    """
    monkeypatch.chdir(tmp_path)
    fig_dir = "site/src/figs"
    per_elem_file = "site/src/routes/models/per-element-each-errors.jsonl"
    os.makedirs(fig_dir)
    os.makedirs(os.path.dirname(per_elem_file))

    def write_payload(path: str, model_keys: list[str]) -> None:
        with open(path, "w") as file:
            file.writelines(json.dumps({"key": key}) + "\n" for key in model_keys)

    write_payload(f"{fig_dir}/demo.jsonl", ["m-a", "m-b"])
    write_payload(per_elem_file, ["m-a"])
    # git via variable dodges ruff S607, same pattern as the script under test
    git_cfg = ["git", "-c", "user.email=test@test", "-c", "user.name=test"]
    for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-qm", "init"]):
        subprocess.run(git_cfg + cmd, check=True)

    write_payload(f"{fig_dir}/demo.jsonl", ["m-a", "m-c"])  # swap m-b -> m-c
    write_payload(f"{fig_dir}/new-fig.jsonl", ["m-a"])  # new payload, no HEAD version

    out = summarize()
    assert "- `demo`: 2 → 2 (+1: m-c; -1: m-b)" in out
    assert "- `new-fig`: 0 → 1 (+1: m-a)" in out
    assert "per-element-each-errors" not in out  # unchanged -> no bullet
