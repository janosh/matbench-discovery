"""Tests for the data-only figure payload exporter in matbench_discovery.figs."""

from __future__ import annotations

import gzip
import json
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from matbench_discovery import figs
from matbench_discovery.cli import cli_args
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
        (None, []),  # None -> [] so a missing trace field stays JSON-serializable
    ],
)
def test_round_list(values: list[Any] | None, expected: list[Any]) -> None:
    """round_list rounds floats, nulls non-finite values, keeps ints/strings."""
    assert figs.round_list(values) == expected


def test_lttb_keeps_endpoints_and_count() -> None:
    """LTTB down-samples to the requested count while preserving first/last points."""
    x = np.linspace(0, 10, 1000)
    y = np.sin(x)
    ds_x, ds_y = figs.lttb(x, y, 50)
    # LTTB targets n_out points but may drop a duplicate bucket index after dedup
    assert 48 <= len(ds_x) <= 50
    assert ds_x[0] == x[0]
    assert ds_x[-1] == x[-1]
    assert ds_y[0] == y[0]
    assert ds_y[-1] == y[-1]


def test_histogram_bins_raw_values() -> None:
    """Histogram returns bin centers, integer counts and bar width, dropping NaNs."""
    values = [0.0, 0.1, 0.1, 0.9, float("nan")]
    result = figs.histogram(values, bins=10, value_range=(0, 1))
    assert set(result) == {"x", "y", "bar_width"}
    assert len(result["x"]) == len(result["y"]) == 10
    assert sum(result["y"]) == 4  # NaN dropped
    assert result["bar_width"] == pytest.approx(0.1)
    assert result["x"][0] == pytest.approx(0.05)  # first bin center


def test_sankey_flow_canonicalization() -> None:
    """Sankey flow data drops unused nodes and sorts links deterministically."""
    # "X" is unreferenced; links are deliberately non-canonical.
    flow_data = {
        "labels": ["A", "B", "X", "C"],
        "source_indices": [1, 0],
        "target_indices": [3, 3],
        "value": [4.0, 3.0],
    }
    assert figs.sankey_payload_from_flow(flow_data) == {
        "labels": ["A", "B", "C"],
        "source": [0, 1],
        "target": [2, 2],
        "value": [3.0, 4.0],
    }


def test_write_json_gz_roundtrip(tmp_path: Path) -> None:
    """write_json_gz writes deterministic gzipped JSON parseable back unchanged."""
    payload = {"models": [{"label": "demo", "x": [1, 2], "y": [3.5, 4.5]}]}
    out_path = f"{tmp_path}/sub/dir/demo.json.gz"  # also creates parent dirs
    size = figs.write_json_gz(out_path, payload)
    assert size > 0
    with gzip.open(out_path) as file:
        assert json.load(file) == payload
    # deterministic output: same payload -> identical bytes (mtime pinned)
    with open(out_path, "rb") as file:
        first_bytes = file.read()
    figs.write_json_gz(out_path, payload)
    with open(out_path, "rb") as file:
        assert file.read() == first_bytes


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
    payload = {"y": [1, 2]}
    path = f"{tmp_path}/demo.json.gz"
    with open(path, "wb") as file:
        file.write(existing_bytes)

    figs.write_json_gz(path, payload)
    with open(path, "rb") as file:
        written_bytes = file.read()
    if preserve_existing:
        assert written_bytes == existing_bytes  # content-equal -> bytes untouched
        return
    with gzip.open(path) as file:
        assert json.load(file) == payload


# === multi-model site payload writing (JSONL: full-run vs subset-run splice) ===
@pytest.fixture
def site_fig_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect write_site_payload's IO to a temp dir and default to a full --models run
    (subset tests narrow cli_args.models themselves).
    """
    import matbench_discovery

    monkeypatch.setattr(matbench_discovery, "SITE_FIG_DATA", str(tmp_path))
    monkeypatch.setattr(cli_args, "models", list(Model.active()))
    return tmp_path


def test_read_jsonl_payload_roundtrip(tmp_path: Path) -> None:
    """read_jsonl_payload routes the lone {"_base": ...} line to shared fields and every
    other line to models, regardless of line order, skipping blank lines.
    """
    path = f"{tmp_path}/fig.jsonl"
    lines = [  # base line need not come first; trailing blank line is ignored
        '{"key":"m-b","y":[2.0]}',
        '{"_base":{"x":[0.0,1.0]}}',
        '{"key":"m-a","y":[1.0]}',
        "",
    ]
    with open(path, "w") as file:
        file.write("\n".join(lines) + "\n")
    restored = figs.read_jsonl_payload(path)
    assert restored["x"] == [0.0, 1.0]
    by_key = {model["key"]: model["y"] for model in restored["models"]}
    assert by_key == {"m-a": [1.0], "m-b": [2.0]}


def test_write_site_payload_full_run_writes_jsonl_and_strips(
    site_fig_dir: Path,
) -> None:
    """Full runs write a lone _base line for shared fields + one line per model, sorted
    by id and stripped of color/visible so lines stay position-independent.
    """
    model_a, model_b = list(Model.active())[:2]
    figs.write_site_payload(
        "demo",
        {
            "shared": [1],
            "models": [
                {"key": model_b.key, "y": [2], "color": "#x", "visible": False},
                {"key": model_a.key, "y": [1]},
            ],
        },
    )
    path = f"{site_fig_dir}/demo.jsonl"
    with open(path) as file:
        raw = file.read().splitlines()
    assert json.loads(raw[0]) == {"_base": {"shared": [1]}}  # _base line written first
    models = [json.loads(line) for line in raw[1:]]
    # presentation stripped + models sorted by id
    assert all("color" not in mdl and "visible" not in mdl for mdl in models)
    assert [mdl["key"] for mdl in models] == sorted([model_a.key, model_b.key])
    restored = figs.read_jsonl_payload(path)
    assert restored["shared"] == [1]
    by_key = {mdl["key"]: mdl["y"] for mdl in restored["models"]}
    assert by_key == {model_a.key: [1], model_b.key: [2]}


def test_write_site_payload_no_base_line_without_shared(site_fig_dir: Path) -> None:
    """No _base line is written when the payload carries no fields beyond models."""
    model_a = next(iter(Model.active()))
    figs.write_site_payload("demo", {"models": [{"key": model_a.key, "y": [1]}]})
    with open(f"{site_fig_dir}/demo.jsonl") as file:
        raw = file.read().splitlines()
    assert all("_base" not in line for line in raw)


def test_write_site_payload_full_run_prunes_dropped_models(site_fig_dir: Path) -> None:
    """A later full run rewrites the whole roster, dropping models no longer present."""
    model_a, model_b = list(Model.active())[:2]
    figs.write_site_payload(
        "demo", {"models": [{"key": model_a.key}, {"key": model_b.key}]}
    )
    figs.write_site_payload("demo", {"models": [{"key": model_b.key}]})  # shrunk roster
    reread = figs.read_jsonl_payload(f"{site_fig_dir}/demo.jsonl")
    assert {model["key"] for model in reread["models"]} == {model_b.key}


@pytest.mark.parametrize("id_field", ["key", "label"])
def test_write_site_payload_subset_run_splices(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch, id_field: str
) -> None:
    """Subset runs (--models) splice fresh entries into the committed file by id_field
    (key- or label-keyed payloads): update their own line, add new ones, leave every
    other model untouched.
    """
    figs.write_site_payload(
        "demo",
        {"models": [{id_field: "m-a", "y": [0]}, {id_field: "m-b", "y": [1]}]},
        id_field=id_field,
    )
    monkeypatch.setattr(cli_args, "models", list(Model.active())[:1])  # subset run
    figs.write_site_payload(
        "demo",
        {"models": [{id_field: "m-a", "y": [9]}, {id_field: "m-c", "y": [3]}]},
        id_field=id_field,
    )
    reread = figs.read_jsonl_payload(f"{site_fig_dir}/demo.jsonl")
    by_id = {model[id_field]: model["y"] for model in reread["models"]}
    assert by_id == {"m-a": [9], "m-b": [1], "m-c": [3]}


@pytest.mark.parametrize("id_field", ["key", "label"])
def test_write_site_payload_subset_run_prunes_inactive_models(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch, id_field: str
) -> None:
    """Subset runs drop committed entries of superseded/inactive models (by key or
    label) while preserving unknown reference lines like 'Test set standard deviation'.
    """
    inactive = next(model for model in Model if not model.is_active)
    stale_id = getattr(inactive, id_field)
    figs.write_site_payload(
        "demo",
        {
            "models": [
                {id_field: "m-a", "y": [0]},
                {id_field: stale_id, "y": [1]},
                {id_field: "Test set standard deviation", "y": [2]},
            ]
        },
        id_field=id_field,
    )
    monkeypatch.setattr(cli_args, "models", list(Model.active())[:1])  # subset run
    figs.write_site_payload(
        "demo", {"models": [{id_field: "m-a", "y": [9]}]}, id_field=id_field
    )
    reread = figs.read_jsonl_payload(f"{site_fig_dir}/demo.jsonl")
    by_id = {model[id_field]: model["y"] for model in reread["models"]}
    assert by_id == {"m-a": [9], "Test set standard deviation": [2]}


def test_write_site_payload_subset_noop_is_byte_identical(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Output is deterministic: a subset run that re-supplies a model's existing data
    leaves the committed file byte-identical (so a no-op refresh opens no churn PR) and
    preserves both the shared _base line and the untouched model.
    """
    path = f"{site_fig_dir}/demo.jsonl"
    payload = {
        "shared": [1],
        "models": [{"key": "m-b", "y": [2]}, {"key": "m-a", "y": [1]}],
    }
    figs.write_site_payload("demo", payload)
    with open(path, "rb") as file:
        full = file.read()
    monkeypatch.setattr(cli_args, "models", list(Model.active())[:1])  # subset run
    figs.write_site_payload(
        "demo", {"shared": [1], "models": [{"key": "m-a", "y": [1]}]}
    )
    with open(path, "rb") as file:
        assert file.read() == full  # m-a unchanged, m-b + _base preserved


def test_write_site_payload_subset_preserves_committed_base(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Subset runs keep the committed shared _base; fresh shared fields are ignored
    (shared data is model-independent - changing it needs a full run).
    """
    figs.write_site_payload(
        "demo", {"shared": [1, 2], "models": [{"key": "m-a", "y": [0]}]}
    )
    monkeypatch.setattr(cli_args, "models", list(Model.active())[:1])  # subset run
    figs.write_site_payload(
        "demo", {"shared": [9, 9], "models": [{"key": "m-a", "y": [1]}]}
    )
    restored = figs.read_jsonl_payload(f"{site_fig_dir}/demo.jsonl")
    assert restored["shared"] == [1, 2]  # committed _base preserved, fresh ignored
    assert restored["models"] == [{"key": "m-a", "y": [1]}]  # model line updated


@pytest.mark.usefixtures("site_fig_dir")
def test_write_site_payload_subset_run_requires_existing_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Subset runs can't create a payload from scratch (would ship a partial roster)."""
    monkeypatch.setattr(cli_args, "models", [next(iter(Model.active()))])
    with pytest.raises(FileNotFoundError, match="splice into an existing"):
        figs.write_site_payload("missing", {"models": []})


def test_payload_writers_reject_nan(site_fig_dir: Path) -> None:
    """Both writers fail loudly on NaN (invalid JSON) instead of writing NaN tokens."""
    nan_payload = {"models": [{"key": "m-a", "y": [float("nan")]}]}
    with pytest.raises(ValueError, match="Out of range float"):
        figs.write_json_gz(f"{site_fig_dir}/bad.json.gz", {"y": [float("nan")]})
    with pytest.raises(ValueError, match="Out of range float"):
        figs.write_site_payload("bad", nan_payload)
