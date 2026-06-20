"""Tests for the update-site-figs PR-body roster summarizer."""

import json

import pytest

from scripts.summarize_payload_changes import format_change, roster


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
