"""Tests for the data-only figure payload exporter in matbench_discovery.figs."""

from __future__ import annotations

import base64
import gzip
import json
from typing import TYPE_CHECKING, Any

import numpy as np
import plotly.graph_objects as go
import pytest

from matbench_discovery import figs
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path


def make_bdata(values: list[float], dtype: str = "f8") -> dict[str, str]:
    """Encode a numeric list as a Plotly base64 typed-array dict."""
    arr = np.asarray(values, dtype=np.dtype(dtype))
    return {"dtype": dtype, "bdata": base64.b64encode(arr.tobytes()).decode()}


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (make_bdata([1.0, 2.0, 3.0]), [1.0, 2.0, 3.0]),
        (make_bdata([4, 5, 6], dtype="i4"), [4, 5, 6]),
        ([7.0, 8.0], [7.0, 8.0]),
        (None, None),
    ],
)
def test_decode_array(
    value: dict[str, str] | list[float] | None, expected: list[float] | None
) -> None:
    """decode_array handles base64 typed arrays, plain lists, and None."""
    result = figs.decode_array(value)
    if expected is None:
        assert result is None
    else:
        assert result is not None
        np.testing.assert_allclose(result, expected, rtol=0, atol=0)


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


def test_trace_helpers_extract_xy_and_color() -> None:
    """trace_xy/trace_color/trace_payload read plotly trace objects."""
    trace = go.Scatter(
        name="demo", x=[1, 2, 3], y=[4.0, 5.0, 6.0], line=dict(color="#123456")
    )
    x, y = figs.trace_xy(trace)
    np.testing.assert_array_equal(x, [1, 2, 3])
    np.testing.assert_array_equal(y, [4.0, 5.0, 6.0])
    assert figs.trace_color(trace) == "#123456"
    assert figs.trace_payload(trace, x=False) == {
        "label": "demo",
        "color": "#123456",
        "y": [4.0, 5.0, 6.0],
    }
    assert figs.trace_payload(trace)["x"] == [1, 2, 3]

    bar = go.Bar(x=[1], y=[2], marker=dict(color="#abcdef"))
    assert figs.trace_color(bar) == "#abcdef"


def test_sunburst_data_extracts_flat_arrays() -> None:
    """sunburst_data returns the flat labels/parents/values/ids arrays unchanged."""
    ids = ["cubic", "cubic/225", "cubic/221", "hexagonal"]
    labels = ["cubic", "225", "221", "hexagonal"]
    parents = ["", "cubic", "cubic", ""]
    values = [10, 6, 4, 5]
    fig = go.Figure(
        go.Sunburst(
            branchvalues="total", ids=ids, labels=labels, parents=parents, values=values
        )
    )
    result = figs.sunburst_data(fig)
    assert result == {
        "labels": labels,
        "parents": parents,
        "values": values,
        "ids": ids,
    }


def test_sankey_data_from_sankey_trace() -> None:
    """sankey_data drops unreferenced nodes, reindexes links onto the kept ones and
    canonicalizes link order (payload bytes must not depend on input link order).
    """
    fig = go.Figure(
        go.Sankey(
            # "X" (index 2) is unreferenced -> dropped; "C" reindexed 3 -> 2.
            # links given in non-canonical order to exercise the link sorting
            node=dict(label=["A", "B", "X", "C"]),
            link=dict(source=[1, 0], target=[3, 3], value=[4.0, 3.0]),
        )
    )
    assert figs.sankey_data(fig) == {
        "labels": ["A", "B", "C"],
        "source": [0, 1],
        "target": [2, 2],
        "value": [3.0, 4.0],
    }


@pytest.mark.parametrize(
    ("converter", "trace_type"),
    [(figs.sunburst_data, "sunburst"), (figs.sankey_data, "sankey")],
)
def test_converters_require_matching_trace(
    converter: Callable[[go.Figure | dict[str, Any]], dict[str, Any]], trace_type: str
) -> None:
    """sunburst_data/sankey_data raise on figures without their trace type."""
    with pytest.raises(ValueError, match=f"no {trace_type} trace"):
        converter(go.Figure(go.Scatter(x=[1], y=[2])))


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


def test_write_json_gz_rejects_nan(tmp_path: Path) -> None:
    """write_json_gz refuses NaN values (invalid JSON) instead of writing them."""
    with pytest.raises(ValueError, match="Out of range float values"):
        figs.write_json_gz(f"{tmp_path}/bad.json.gz", {"y": [float("nan")]})


# === multi-model site payload writing (JSONL: full-run vs subset-run splice) ===
@pytest.fixture
def site_fig_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect write_site_payload's IO to a temp dir."""
    import matbench_discovery

    monkeypatch.setattr(matbench_discovery, "SITE_FIG_DATA", str(tmp_path))
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
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Full runs write a lone _base line for shared fields + one line per model, sorted
    by id and stripped of color/visible so lines stay position-independent.
    """
    model_a, model_b = list(Model.active())[:2]
    monkeypatch.setattr(cli_args, "models", list(Model.active()))  # full run
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


def test_write_site_payload_no_base_line_without_shared(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No _base line is written when the payload carries no fields beyond models."""
    model_a = next(iter(Model.active()))
    monkeypatch.setattr(cli_args, "models", list(Model.active()))
    figs.write_site_payload("demo", {"models": [{"key": model_a.key, "y": [1]}]})
    with open(f"{site_fig_dir}/demo.jsonl") as file:
        raw = file.read().splitlines()
    assert all("_base" not in line for line in raw)


def test_write_site_payload_full_run_prunes_dropped_models(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A later full run rewrites the whole roster, dropping models no longer present."""
    model_a, model_b = list(Model.active())[:2]
    monkeypatch.setattr(cli_args, "models", list(Model.active()))
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
    monkeypatch.setattr(cli_args, "models", list(Model.active()))  # full run first
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
    monkeypatch.setattr(cli_args, "models", list(Model.active()))  # full run first
    figs.write_site_payload("demo", payload)
    with open(path, "rb") as file:
        full = file.read()
    monkeypatch.setattr(cli_args, "models", list(Model.active())[:1])  # subset run
    figs.write_site_payload(
        "demo", {"shared": [1], "models": [{"key": "m-a", "y": [1]}]}
    )
    with open(path, "rb") as file:
        assert file.read() == full  # m-a unchanged, m-b + _base preserved


@pytest.mark.usefixtures("site_fig_dir")
def test_write_site_payload_subset_run_requires_existing_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Subset runs can't create a payload from scratch (would ship a partial roster)."""
    monkeypatch.setattr(cli_args, "models", [next(iter(Model.active()))])
    with pytest.raises(FileNotFoundError, match="splice into an existing"):
        figs.write_site_payload("missing", {"models": []})


@pytest.mark.usefixtures("site_fig_dir")
def test_write_site_payload_rejects_nan(monkeypatch: pytest.MonkeyPatch) -> None:
    """NaN values (invalid JSON) fail loudly instead of writing literal NaN tokens."""
    model_a = next(iter(Model.active()))
    monkeypatch.setattr(cli_args, "models", list(Model.active()))
    with pytest.raises(ValueError, match="Out of range float"):
        figs.write_site_payload(
            "bad", {"models": [{"key": model_a.key, "y": [float("nan")]}]}
        )
