"""Tests for deterministic data-only figure payload export."""

from __future__ import annotations

import gzip
import json
from typing import TYPE_CHECKING, Any
from unittest.mock import patch

import numpy as np
import pytest

from matbench_discovery import figs, payload_numerics

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


def test_discovery_provenance_hashes_only_caller_recipe_sources() -> None:
    """Discovery identity does not hash data.py or invent unused recipe sources."""
    with patch.object(figs, "build_payload_provenance", return_value={}) as builder:
        figs.build_discovery_payload_provenance(
            generator=__file__, test_subset="test", parameters={}
        )
        assert builder.call_args.kwargs["source_files"] == {}
        figs.build_discovery_payload_provenance(
            generator=__file__,
            test_subset="test",
            parameters={},
            source_files={"payload_numerics": payload_numerics.__file__},
        )
    assert builder.call_args.kwargs["source_files"] == {
        "payload_numerics": payload_numerics.__file__
    }


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


def make_metadata(
    tmp_path: Path,
    *,
    parameters: dict[str, object] | None = None,
) -> dict[str, Any]:
    """Build valid test identity and audit metadata from exact source bytes."""
    benchmark_path = tmp_path / "benchmark.csv"
    benchmark_path.write_text("benchmark")
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
    input_content: str | None = None,
) -> dict[str, Any]:
    """Build one valid model record with a path-free input manifest."""
    content = input_content or f"input-{model_key}"
    input_path = tmp_path / f"{model_key}.csv"
    input_path.write_text(content)
    return {
        "model_key": model_key,
        "label": model_key.upper(),
        "input_artifacts": [figs.artifact_manifest("predictions", str(input_path))],
        "value": value,
    }


def write_test_payload(
    path: str,
    tmp_path: Path,
    models: list[dict[str, Any]],
    *,
    shared: object = 1,
    metadata: dict[str, Any] | None = None,
    target_keys: set[str] | None = None,
) -> dict[str, Any]:
    """Write and read one complete test payload."""
    figs.write_jsonl_payload(
        path,
        {
            **(metadata or make_metadata(tmp_path)),
            "shared": shared,
            "models": models,
        },
        target_keys=target_keys,
    )
    return figs.read_jsonl_payload(path)


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
    assert set(json.loads(lines[0])["_base"]) == {
        "schema_version",
        "identity",
        "audit",
        "derived",
    }
    assert [json.loads(line)["model_key"] for line in lines[1:]] == [
        "model-a",
        "model-b",
    ]
    assert all(line == figs.canonical_json(json.loads(line)) for line in lines)
    assert all("color" not in line and "visible" not in line for line in lines[1:])
    write_test_payload(path, tmp_path, models)
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


def test_writer_replaces_targeted_records_or_full_payload(tmp_path: Path) -> None:
    """Targeted writes preserve peers while full writes replace the payload."""
    path = f"{tmp_path}/payload.jsonl"
    model_a = make_model(tmp_path, "model-a", [1])
    model_b = make_model(tmp_path, "model-b", [2])
    before = write_test_payload(path, tmp_path, [model_a, model_b])
    changed_b = make_model(
        tmp_path, "model-b", [9], input_content="changed-model-b-input"
    )
    with pytest.raises(ValueError, match="unselected model keys"):
        write_test_payload(
            path, tmp_path, [model_a, changed_b], target_keys={"model-a"}
        )
    changed_a = make_model(
        tmp_path, "model-a", [9], input_content="changed-model-a-input"
    )
    metadata = make_metadata(tmp_path)
    metadata["audit"]["runtime"]["python"] = "0.0.0"
    after = write_test_payload(
        path,
        tmp_path,
        [changed_a],
        metadata=metadata,
        target_keys={"model-a"},
    )
    assert after["models"] == [changed_a, model_b]
    assert before | {"audit": metadata["audit"], "models": after["models"]} == after
    with pytest.raises(ValueError, match="identity changed"):
        write_test_payload(
            path,
            tmp_path,
            [changed_a],
            metadata=make_metadata(tmp_path, parameters={"changed": True}),
            target_keys={"model-a"},
        )
    with pytest.raises(ValueError, match="shared data changed"):
        write_test_payload(
            path, tmp_path, [changed_a], shared=999, target_keys={"model-a"}
        )
    removed = write_test_payload(path, tmp_path, [], target_keys={"model-a"})
    assert removed["models"] == [model_b]
    model_c = make_model(tmp_path, "model-c", [3])
    metadata = make_metadata(tmp_path, parameters={"version": 2})
    with pytest.raises(ValueError, match="identity and roster changed together"):
        write_test_payload(
            path, tmp_path, [model_b, model_c], shared=2, metadata=metadata
        )
    roster_updated = write_test_payload(path, tmp_path, [model_b, model_c])
    assert roster_updated["models"] == [model_b, model_c]
    replaced = write_test_payload(
        path, tmp_path, [model_b, model_c], shared=2, metadata=metadata
    )
    assert replaced["models"] == [model_b, model_c]
    assert replaced["identity"] == metadata["identity"]
    assert replaced["shared"] == 2
    metadata["audit"]["source_commit"] = "0" * 40
    with pytest.raises(ValueError, match="audit must contain runtime"):
        write_test_payload(
            path, tmp_path, [model_b, model_c], shared=2, metadata=metadata
        )


@pytest.mark.parametrize("target_keys", [None, {"model-a"}], ids=["full", "targeted"])
def test_writer_allows_labels_but_rejects_data_changes_with_unchanged_inputs(
    tmp_path: Path, target_keys: set[str] | None
) -> None:
    """Stable computation inputs permit labels but not derived-data changes."""
    path = f"{tmp_path}/payload.jsonl"
    model = make_model(tmp_path, "model-a", [1])
    write_test_payload(path, tmp_path, [model])
    relabeled = write_test_payload(
        path, tmp_path, [model | {"label": "Renamed"}], target_keys=target_keys
    )
    assert relabeled["models"][0]["label"] == "Renamed"
    with pytest.raises(ValueError, match="changed data with unchanged inputs"):
        write_test_payload(
            path,
            tmp_path,
            [model | {"value": [2]}],
            target_keys=target_keys,
        )


def test_payload_writers_reject_nan(tmp_path: Path) -> None:
    """Both payload writers reject non-standard NaN JSON tokens."""
    with pytest.raises(ValueError, match="Out of range float"):
        figs.write_json_gz(f"{tmp_path}/bad.json.gz", {"y": [float("nan")]})
    model = make_model(tmp_path, "model-a", [float("nan")])
    with pytest.raises(ValueError, match="Out of range float"):
        write_test_payload(f"{tmp_path}/bad.jsonl", tmp_path, [model])
