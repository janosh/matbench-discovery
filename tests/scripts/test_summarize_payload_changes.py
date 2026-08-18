"""Tests for structural and numerical payload-change reporting."""

import json
import os
import subprocess
from typing import TYPE_CHECKING, Any

import pytest

from scripts.summarize_payload_changes import _parse_payload, compare_values, summarize

if TYPE_CHECKING:
    from pathlib import Path


def test_compare_values_reports_numeric_and_structural_deltas() -> None:
    """Matching numeric leaves get exact errors; shape/type changes stay separate."""
    delta = compare_values(
        {"values": [0, 2], "flag": True, "missing": [1]},
        {"values": [0, 4], "flag": False, "added": "x"},
    )
    assert delta.numeric_leaf_count == 2
    assert delta.max_absolute_error == 2
    assert delta.max_relative_error == 0.5
    assert delta.structural_paths == [
        "$.added (added)",
        "$.flag",
        "$.missing (removed)",
    ]

    identity_old = {
        "model_key": "old",
        "input_artifacts": [{"role": "x", "sha256": "0" * 64, "size": 1}],
        "value": 2,
    }
    identity_new = {
        "model_key": "new",
        "input_artifacts": [{"role": "x", "sha256": "1" * 64, "size": 99}],
        "value": 3,
    }
    identity_delta = compare_values(identity_old, identity_new)
    assert identity_delta.numeric_leaf_count == 1
    assert identity_delta.max_absolute_error == 1
    assert identity_delta.max_relative_error == 1 / 3
    assert identity_delta.structural_paths == [
        "$.input_artifacts[0].sha256",
        "$.input_artifacts[0].size",
        "$.model_key",
    ]
    shape_delta = compare_values([1, 2], [1, 99, 3])
    assert shape_delta.structural_paths == ["$ (length 2 -> 3)"]
    assert shape_delta.numeric_leaf_count == 2
    assert shape_delta.max_absolute_error == 97
    assert shape_delta.max_relative_error == 97 / 99


def test_parse_payload_maps_historical_identity_fields() -> None:
    """The reporter maps historical key and label records to current model keys."""
    _base, models = _parse_payload(
        '{"key":"legacy-key"}\n{"label":"Old display name"}\n',
        {"Old display name": "current-key"},
    )
    assert set(models) == {"current-key", "legacy-key"}


def payload_text(
    models: dict[str, int],
    *,
    shared: int = 1,
    labels: dict[str, str] | None = None,
    identity: dict[str, object] | None = None,
    audit: dict[str, object] | None = None,
) -> str:
    """Serialize a minimal schema-v2-like JSONL fixture for report tests."""
    base = {
        "schema_version": 2,
        "identity": identity or {"recipe": {"sha256": "a" * 64}},
        "audit": audit or {"runtime": {"python": "3.14", "packages": {}}},
        "derived": {"shared": shared},
    }
    records = [{"_base": base}] + [
        {
            "model_key": model_key,
            "label": (labels or {}).get(model_key, model_key),
            "value": value,
        }
        for model_key, value in models.items()
    ]
    return "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)


def write_payload(
    path: str,
    models: dict[str, int],
    *,
    shared: int = 1,
    labels: dict[str, str] | None = None,
    identity: dict[str, object] | None = None,
    audit: dict[str, object] | None = None,
) -> None:
    """Write one minimal JSONL payload fixture."""
    with open(path, "w", encoding="utf-8") as file:
        file.write(
            payload_text(
                models,
                shared=shared,
                labels=labels,
                identity=identity,
                audit=audit,
            )
        )


def test_summarize_diffs_working_tree_against_head(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end report includes record, roster, identity and numerical deltas."""
    monkeypatch.chdir(tmp_path)
    fig_dir = "site/src/figs"
    per_element_file = "site/src/routes/models/per-element-each-errors.jsonl"
    os.makedirs(fig_dir)
    os.makedirs(os.path.dirname(per_element_file))
    demo_path = f"{fig_dir}/demo.jsonl"
    deleted_path = f"{fig_dir}/deleted.jsonl"
    identity: dict[str, Any] = {
        "benchmark_inputs": [{"role": "bench", "sha256": "0" * 64, "size": 1}],
        "parameters": {"rounding": 3},
        "recipe": {"sha256": "a" * 64},
    }
    audit: dict[str, Any] = {"runtime": {"python": "3.14", "packages": {}}}
    write_payload(
        demo_path,
        {"model-a": 1, "model-b": 2},
        identity=identity,
        audit=audit,
    )
    write_payload(per_element_file, {"model-a": 1})
    write_payload(deleted_path, {"model-a": 1})

    git_config = ["git", "-c", "user.email=test@test", "-c", "user.name=test"]
    for command in (["init", "-q"], ["add", "-A"], ["commit", "-qm", "init"]):
        subprocess.run(git_config + command, check=True)

    changed_identity = json.loads(json.dumps(identity))
    changed_identity["benchmark_inputs"][0]["sha256"] = "1" * 64
    changed_identity["parameters"]["rounding"] = 4
    changed_audit = json.loads(json.dumps(audit))
    changed_audit["runtime"]["python"] = "3.15"
    write_payload(
        demo_path,
        {"model-a": 1, "model-c": 4},
        shared=3,
        labels={"model-a": "Renamed A"},
        identity=changed_identity,
        audit=changed_audit,
    )
    write_payload(f"{fig_dir}/new-fig.jsonl", {"model-a": 1})
    write_payload(per_element_file, {"model-a": 2})
    os.remove(deleted_path)

    report = summarize()
    assert "| `demo` | `<shared>` |" in report
    assert "| `demo` | `model-a` | `$.label` | 1 | 0 | 0 |" in report
    assert "`demo` roster 2 -> 2; added: model-c; removed: model-b" in report
    assert "`new-fig` roster 0 -> 1; added: model-a; removed: —" in report
    assert "`deleted` roster 1 -> 0; added: —; removed: model-a" in report
    assert "Max symmetric relative error" in report
    assert "benchmark_inputs, parameters" in report
    assert "| `per-element-each-errors` | `model-a` |" in report
