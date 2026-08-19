"""Tests for key-addressed payload-change reporting."""

import json
import os
import subprocess
from typing import TYPE_CHECKING

import pytest

from scripts.summarize_payload_changes import _parse_payload, summarize

if TYPE_CHECKING:
    from pathlib import Path


def test_parse_payload_requires_unique_model_keys() -> None:
    """The reporter requires one immutable key per model record."""
    with pytest.raises(TypeError, match="no model_key"):
        _parse_payload('{"label":"Display name"}\n')
    record = '{"model_key":"duplicate","label":"Duplicate"}\n'
    with pytest.raises(ValueError, match="Duplicate payload model_key"):
        _parse_payload(record * 2)


def write_payload(
    path: str,
    models: dict[str, int],
    *,
    shared: int = 1,
    identity_version: float = 1,
) -> None:
    """Write one minimal JSONL payload fixture."""
    records = [
        {
            "_base": {
                "schema_version": 2,
                "identity": {"parameters": {"version": identity_version}},
                "derived": {"shared": shared},
            }
        },
        *(
            {"model_key": model_key, "label": model_key, "value": value}
            for model_key, value in models.items()
        ),
    ]
    with open(path, "w", encoding="utf-8") as file:
        file.write(
            "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)
        )


def test_summarize_diffs_working_tree_against_head(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The report identifies shared, added, removed, and updated records."""
    monkeypatch.chdir(tmp_path)
    fig_dir = "site/src/figs"
    route_dir = "site/src/routes/models"
    os.makedirs(fig_dir)
    os.makedirs(route_dir)
    demo_path = f"{fig_dir}/demo.jsonl"
    deleted_path = f"{fig_dir}/deleted.jsonl"
    route_path = f"{route_dir}/per-element-each-errors.jsonl"
    write_payload(demo_path, {"model-a": 1, "model-b": 2})
    write_payload(deleted_path, {"model-a": 1})
    write_payload(route_path, {"model-a": 1})

    git_config = ["git", "-c", "user.email=test@test", "-c", "user.name=test"]
    for command in (["init", "-q"], ["add", "-A"], ["commit", "-qm", "init"]):
        subprocess.run(git_config + command, check=True)

    write_payload(
        demo_path,
        {"model-a": 3, "model-c": 4},
        shared=2,
        identity_version=1.0,
    )
    write_payload(f"{fig_dir}/new.jsonl", {"model-a": 1})
    write_payload(route_path, {"model-a": 2})
    os.remove(deleted_path)

    report = summarize()
    assert report == (
        "### Payload changes\n\n"
        "- `deleted`: identity; shared data; removed: model-a\n"
        "- `demo`: identity (parameters); shared data (max |Δ|=1, max rel=0.5); "
        "added: model-c; removed: model-b; updated: model-a "
        "(max |Δ|=2, max rel=0.666667)\n"
        "- `new`: identity; shared data; added: model-a\n"
        "- `per-element-each-errors`: updated: model-a "
        "(max |Δ|=1, max rel=0.5)"
    )
