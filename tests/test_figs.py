"""Tests for deterministic data-only figure payload export."""

from __future__ import annotations

import gzip
import json
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any
from unittest.mock import patch

import numpy as np
import pytest

from matbench_discovery import figs, payload_numerics
from matbench_discovery.enums import Model

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([1, 2, 3], [1, 2, 3]),
        ([1.23456789, 2.0], [1.23457, 2.0]),
        ([1.0, float("nan"), float("inf")], [1.0, None, None]),
        (["Fe", "Co"], ["Fe", "Co"]),
        (None, []),
    ],
)
def test_round_list(values: list[Any] | None, expected: list[Any]) -> None:
    """round_list rounds floats, nulls non-finite values, and keeps other scalars."""
    assert payload_numerics.round_list(values) == expected


@pytest.mark.parametrize(
    ("n_values", "n_out", "expected"),
    [(6, 4, [0, 1, 3, 5]), (3, 5, [0, 1, 2]), (1, 5, [0])],
)
def test_evenly_spaced_indices_are_distinct_and_bounded(
    n_values: int, n_out: int, expected: list[int]
) -> None:
    """Index sampling caps output at the number of distinct available entries."""
    assert payload_numerics.evenly_spaced_indices(n_values, n_out).tolist() == expected


@pytest.mark.parametrize(
    ("n_values", "n_out"), [(6, 5), (10, 5), (100, 20), (1000, 50)]
)
def test_lttb_keeps_endpoints_and_exact_count(n_values: int, n_out: int) -> None:
    """LTTB keeps both endpoints and returns exactly the requested sample count."""
    x_values = np.linspace(0, 10, n_values)
    y_values = np.sin(x_values)
    downsampled_x, downsampled_y = payload_numerics.lttb(x_values, y_values, n_out)
    assert len(downsampled_x) == len(downsampled_y) == n_out
    assert np.all(np.diff(downsampled_x) > 0)
    assert downsampled_x[0] == x_values[0]
    assert downsampled_x[-1] == x_values[-1]
    assert downsampled_y[0] == y_values[0]
    assert downsampled_y[-1] == y_values[-1]


def test_histogram_bins_raw_values() -> None:
    """Histogram returns centers, integer counts and width while dropping NaNs."""
    result = payload_numerics.histogram(
        [0.0, 0.1, 0.1, 0.9, float("nan")], bins=10, value_range=(0, 1)
    )
    assert set(result) == {"x", "y", "bar_width"}
    assert len(result["x"]) == len(result["y"]) == 10
    assert sum(result["y"]) == 4
    assert result["bar_width"] == pytest.approx(0.1)
    assert result["x"][0] == pytest.approx(0.05)


def test_sankey_flow_canonicalization() -> None:
    """Sankey flow data drops unused nodes and sorts links deterministically."""
    flow_data = {
        "labels": ["A", "B", "X", "C"],
        "source_indices": [1, 0],
        "target_indices": [3, 3],
        "value": [4.0, 3.0],
    }
    assert payload_numerics.sankey_payload_from_flow(flow_data) == {
        "labels": ["A", "B", "C"],
        "source": [0, 1],
        "target": [2, 2],
        "value": [3.0, 4.0],
    }


def test_write_json_gz_roundtrip(tmp_path: Path) -> None:
    """write_json_gz writes deterministic gzip bytes parseable unchanged."""
    payload = {"models": [{"label": "demo", "x": [1, 2], "y": [3.5, 4.5]}]}
    out_path = f"{tmp_path}/sub/dir/demo.json.gz"
    assert figs.write_json_gz(out_path, payload) > 0
    with gzip.open(out_path) as file:
        assert json.load(file) == payload
    with open(out_path, "rb") as file:
        first_bytes = file.read()
    figs.write_json_gz(out_path, payload)
    with open(out_path, "rb") as file:
        assert file.read() == first_bytes


def test_artifact_manifest_hashes_and_sizes_one_open_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Artifact sizes come from the file descriptor used for hashing."""
    path = tmp_path / "artifact.csv"
    path.write_text("data")
    monkeypatch.setattr(
        figs.os.path,
        "getsize",
        lambda _path: pytest.fail("artifact manifest must not reopen the path"),
    )
    with patch("builtins.open", wraps=open) as mock_open:
        assert figs.artifact_manifest("input", str(path))["size"] == 4
    mock_open.assert_called_once_with(str(path), "rb")


@pytest.mark.parametrize(
    ("existing_bytes", "preserve_existing"),
    [
        pytest.param(
            gzip.compress(
                json.dumps({"y": [1, 2]}, separators=(",", ":")).encode(),
                compresslevel=1,
            ),
            True,
            id="content-equal-gzip",
        ),
        pytest.param(b"\x1f\x8b\x08\x00", False, id="truncated-gzip"),
        pytest.param(b"this is not a gzip stream", False, id="not-gzip"),
    ],
)
def test_write_json_gz_handles_existing_file(
    tmp_path: Path, existing_bytes: bytes, preserve_existing: bool
) -> None:
    """write_json_gz preserves equivalent gzip bytes and rewrites corrupt files."""
    path = f"{tmp_path}/demo.json.gz"
    with open(path, "wb") as file:
        file.write(existing_bytes)
    figs.write_json_gz(path, {"y": [1, 2]})
    with open(path, "rb") as file:
        written_bytes = file.read()
    if preserve_existing:
        assert written_bytes == existing_bytes
    else:
        with gzip.open(path) as file:
            assert json.load(file) == {"y": [1, 2]}


def make_provenance(
    tmp_path: Path,
    *,
    benchmark_content: str = "benchmark",
    parameters: dict[str, object] | None = None,
) -> dict[str, Any]:
    """Build valid test provenance from exact benchmark and source bytes."""
    benchmark_path = tmp_path / "benchmark.csv"
    benchmark_path.write_text(benchmark_content)
    return figs.build_payload_provenance(
        generator=__file__,
        benchmark_inputs={"benchmark": str(benchmark_path)},
        source_files={},
        parameters=parameters or {"rounding_decimals": 5},
        packages=(),
    )


def make_model(
    tmp_path: Path,
    model_key: str,
    value: object,
    *,
    label: str | None = None,
    input_content: str | None = None,
) -> dict[str, Any]:
    """Build one valid model record with a path-free input manifest."""
    content = input_content or f"input-{model_key}"
    input_path = tmp_path / f"{model_key}.csv"
    input_path.write_text(content)
    return {
        "model_key": model_key,
        "label": label or model_key.upper(),
        "input_artifacts": [figs.artifact_manifest("predictions", str(input_path))],
        "value": value,
    }


def write_test_payload(
    path: str,
    tmp_path: Path,
    models: list[dict[str, Any]],
    *,
    mode: figs.PayloadMode = figs.PayloadMode.migrate_provenance,
    shared: object = 1,
    provenance: dict[str, Any] | None = None,
    key_migration: tuple[str, str] | None = None,
) -> dict[str, Any]:
    """Write and read one complete test payload."""
    figs.write_jsonl_payload(
        path,
        {
            "provenance": provenance or make_provenance(tmp_path),
            "shared": shared,
            "models": models,
        },
        mode=mode,
        key_migration=key_migration,
    )
    return figs.read_jsonl_payload(path)


def test_fingerprint_excludes_source_commit(tmp_path: Path) -> None:
    """Informational source_commit never changes complete computation identity."""
    provenance = make_provenance(tmp_path)
    model = make_model(tmp_path, "model-a", [1])
    base_a = {"provenance": provenance | {"source_commit": "a" * 40}}
    base_b = {"provenance": provenance | {"source_commit": "b" * 40}}
    assert figs.computation_fingerprint(base_a, model) == figs.computation_fingerprint(
        base_b, model
    )


def test_write_jsonl_is_canonical_and_byte_deterministic(tmp_path: Path) -> None:
    """Schema-v2 JSONL is base-first, sorted, canonical, terminated and stable."""
    path = f"{tmp_path}/payload.jsonl"
    models = [
        make_model(tmp_path, "model-b", [2]) | {"color": "red"},
        make_model(tmp_path, "model-a", [1]) | {"visible": False},
    ]
    write_test_payload(path, tmp_path, models)
    with open(path, "rb") as file:
        first_bytes = file.read()
    assert first_bytes.endswith(b"\n")
    assert b"\r" not in first_bytes
    lines = first_bytes.decode().splitlines()
    assert list(json.loads(lines[0])) == ["_base"]
    assert [json.loads(line)["model_key"] for line in lines[1:]] == [
        "model-a",
        "model-b",
    ]
    assert all(line == figs.canonical_json(json.loads(line)) for line in lines)
    assert all("color" not in line and "visible" not in line for line in lines[1:])
    write_test_payload(path, tmp_path, models, mode=figs.PayloadMode.full_roster)
    with open(path, "rb") as file:
        assert file.read() == first_bytes


@pytest.mark.parametrize(
    ("overrides", "count", "match"),
    [
        ({"input_artifacts": []}, 1, "non-empty"),
        ({}, 2, "unique"),
        ({"key": "model-a"}, 1, "Legacy"),
    ],
)
def test_writer_rejects_invalid_model_identity(
    tmp_path: Path, overrides: dict[str, object], count: int, match: str
) -> None:
    """Every model line needs one unique canonical key and non-empty artifacts."""
    model = make_model(tmp_path, "model-a", [1]) | overrides
    with pytest.raises(ValueError, match=match):
        write_test_payload(f"{tmp_path}/bad.jsonl", tmp_path, [model] * count)


def test_label_only_targeted_update_preserves_every_other_byte(tmp_path: Path) -> None:
    """An unchanged fingerprint permits exactly a presentation-label replacement."""
    path = f"{tmp_path}/payload.jsonl"
    old = make_model(tmp_path, "model-a", [1], label="Old label")
    old_record = write_test_payload(path, tmp_path, [old])["models"][0]
    new = make_model(tmp_path, "model-a", [1], label="New label")
    new_record = write_test_payload(
        path, tmp_path, [new], mode=figs.PayloadMode.targeted
    )["models"][0]
    assert old_record | {"label": "New label"} == new_record


def test_unchanged_fingerprint_rejects_derived_mutation(tmp_path: Path) -> None:
    """Fixed computation identity rejects mutations in ordinary and migration runs."""
    path = f"{tmp_path}/payload.jsonl"
    model = make_model(tmp_path, "model-a", [1])
    write_test_payload(path, tmp_path, [model])
    with pytest.raises(ValueError, match="unchanged complete fingerprint"):
        write_test_payload(
            path,
            tmp_path,
            [model | {"value": [2]}],
            mode=figs.PayloadMode.targeted,
        )
    with pytest.raises(ValueError, match="unchanged complete fingerprint"):
        write_test_payload(path, tmp_path, [model | {"value": [999]}])
    with pytest.raises(ValueError, match="shared derived data"):
        write_test_payload(path, tmp_path, [model], shared=999)


@pytest.mark.parametrize(
    "models",
    [
        pytest.param([], id="remove"),
        pytest.param(["model-a", "model-b"], id="add"),
        pytest.param(["model-b"], id="replace"),
    ],
)
def test_provenance_migration_rejects_roster_changes(
    tmp_path: Path, models: list[str]
) -> None:
    """Provenance migration cannot also add, remove, or replace model keys."""
    path = f"{tmp_path}/payload.jsonl"
    write_test_payload(path, tmp_path, [make_model(tmp_path, "model-a", [1])])
    changed_provenance = make_provenance(tmp_path, parameters={"version": 2})
    with pytest.raises(ValueError, match="cannot change the model roster"):
        write_test_payload(
            path,
            tmp_path,
            [make_model(tmp_path, model_key, [1]) for model_key in models],
            provenance=changed_provenance,
        )


def test_targeted_input_change_is_isolated_to_one_model(tmp_path: Path) -> None:
    """A changed model artifact updates only its keyed line and preserves peers/base."""
    path = f"{tmp_path}/payload.jsonl"
    model_a = make_model(tmp_path, "model-a", [1])
    model_b = make_model(tmp_path, "model-b", [2])
    before = write_test_payload(path, tmp_path, [model_a, model_b])
    changed_a = make_model(
        tmp_path, "model-a", [9], input_content="changed-model-a-input"
    )
    after = write_test_payload(
        path, tmp_path, [changed_a], mode=figs.PayloadMode.targeted
    )
    before_by_key = {record["model_key"]: record for record in before["models"]}
    after_by_key = {record["model_key"]: record for record in after["models"]}
    assert after_by_key["model-a"] == changed_a
    assert after_by_key["model-b"] == before_by_key["model-b"]
    assert before | {"models": after["models"]} == after


def test_full_roster_adds_and_removes_records(tmp_path: Path) -> None:
    """Full-roster mode owns lifecycle additions/removals under fixed provenance."""
    path = f"{tmp_path}/payload.jsonl"
    model_a = make_model(tmp_path, "model-a", [1])
    model_b = make_model(tmp_path, "model-b", [2])
    write_test_payload(path, tmp_path, [model_a, model_b])
    model_c = make_model(tmp_path, "model-c", [3])
    updated = write_test_payload(
        path, tmp_path, [model_b, model_c], mode=figs.PayloadMode.full_roster
    )
    assert {record["model_key"] for record in updated["models"]} == {
        "model-b",
        "model-c",
    }

    old = make_model(tmp_path, "old-key", [4], input_content="same-input")
    new = make_model(tmp_path, "new-key", [4], input_content="same-input")
    write_test_payload(path, tmp_path, [old], mode=figs.PayloadMode.full_roster)
    assert write_test_payload(path, tmp_path, [new], mode=figs.PayloadMode.full_roster)[
        "models"
    ] == [new]


@pytest.mark.parametrize("identity_field", ["benchmark", "parameters", "runtime"])
def test_shared_identity_change_requires_provenance_migration(
    tmp_path: Path, identity_field: str
) -> None:
    """Benchmark, parameter, and runtime changes fail outside migration mode."""
    path = f"{tmp_path}/payload.jsonl"
    model = make_model(tmp_path, "model-a", [1])
    old_provenance = make_provenance(tmp_path)
    write_test_payload(path, tmp_path, [model], provenance=old_provenance)
    if identity_field == "benchmark":
        new_provenance = make_provenance(tmp_path, benchmark_content="changed")
    else:
        new_provenance = json.loads(figs.canonical_json(old_provenance))
        if identity_field == "parameters":
            new_provenance["parameters"]["rounding_decimals"] = 4
        else:
            new_provenance["runtime"]["python"] = "0.0.0"
    changed_payload = {
        "provenance": new_provenance,
        "shared": 1,
        "models": [model],
    }
    with pytest.raises(ValueError, match="--migrate-provenance"):
        figs.write_jsonl_payload(
            path, changed_payload, mode=figs.PayloadMode.full_roster
        )
    written_provenance = write_test_payload(
        path, tmp_path, [model], provenance=new_provenance
    )["provenance"]
    assert written_provenance == new_provenance


def test_explicit_key_migration_requires_alias_and_exact_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Key migration accepts only an aliased exact identity rewrite."""
    monkeypatch.setattr(
        Model,
        "from_ref",
        staticmethod(
            lambda model_key: SimpleNamespace(key=model_key, key_aliases=("old-key",))
        ),
    )
    path = f"{tmp_path}/payload.jsonl"
    old = make_model(tmp_path, "old-key", [1], input_content="same-input")
    write_test_payload(path, tmp_path, [old])

    new = make_model(tmp_path, "new-key", [1], input_content="same-input")
    for _attempt in range(2):
        migrated = write_test_payload(
            path,
            tmp_path,
            [new],
            mode=figs.PayloadMode.migrate_model_key,
            key_migration=("old-key", "new-key"),
        )
        assert migrated["models"] == [new]

    path = f"{tmp_path}/peers.jsonl"
    peer = make_model(tmp_path, "peer", [2])
    write_test_payload(path, tmp_path, [old, peer])
    with pytest.raises(ValueError, match="peer record"):
        write_test_payload(
            path,
            tmp_path,
            [new, peer | {"value": [999]}],
            mode=figs.PayloadMode.migrate_model_key,
            key_migration=("old-key", "new-key"),
        )


def test_payload_writers_reject_nan(tmp_path: Path) -> None:
    """Both payload writers reject non-standard NaN JSON tokens."""
    with pytest.raises(ValueError, match="Out of range float"):
        figs.write_json_gz(f"{tmp_path}/bad.json.gz", {"y": [float("nan")]})
    model = make_model(tmp_path, "model-a", [float("nan")])
    with pytest.raises(ValueError, match="Out of range float"):
        write_test_payload(f"{tmp_path}/bad.jsonl", tmp_path, [model])
