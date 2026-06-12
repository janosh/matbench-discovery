"""Tests for the data-only figure payload exporter in matbench_discovery.figs."""

from __future__ import annotations

import base64
import gzip
import json
from typing import TYPE_CHECKING, Any

import numpy as np
import plotly.graph_objects as go
import pytest
from plotly.express.colors import qualitative

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
    """sankey_data drops unreferenced nodes and reindexes links onto the kept ones."""
    fig = go.Figure(
        go.Sankey(
            # "X" (index 2) is unreferenced -> dropped; "C" reindexed 3 -> 2
            node=dict(label=["A", "B", "X", "C"]),
            link=dict(source=[0, 1], target=[3, 3], value=[3.0, 4.0]),
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


# === multi-model site payload writing (full-run overwrite vs subset-run merge) ===
def load_payload(path: str) -> dict[str, Any]:
    """Read back a gzipped JSON payload."""
    with gzip.open(path) as file:
        return json.load(file)


@pytest.fixture
def site_fig_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect write_site_payload's IO to a temp dir."""
    import matbench_discovery

    monkeypatch.setattr(matbench_discovery, "SITE_FIG_DATA", str(tmp_path))
    return tmp_path


def test_write_site_payload_full_run_overwrites_and_styles(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Full runs overwrite wholesale, prune retired models and assign deterministic
    order/colors/visibility.
    """
    model_a, model_b = list(Model.active())[:2]
    # pre-existing payload whose contents must be fully replaced
    figs.write_json_gz(
        f"{site_fig_dir}/demo.json.gz", {"models": [{"key": "stale", "y": [0]}]}
    )
    monkeypatch.setattr(cli_args, "models", list(Model.active()))
    fresh = {
        "shared": "ref-data",
        "models": [
            {"key": model_b.key, "mae": 2.0, "color": "#000", "visible": False},
            {"key": model_a.key, "mae": 1.0, "color": "#fff"},
            {"key": "retired-model", "mae": 0.0},
        ],
    }
    figs.write_site_payload(
        "demo", fresh, sort_key=lambda entry: entry["mae"], assign_colors=True,
        visible_top_n=1,
    )  # fmt: skip
    written = load_payload(f"{site_fig_dir}/demo.json.gz")
    assert written["shared"] == "ref-data"
    # retired model pruned, remaining entries sorted by mae
    assert [entry["key"] for entry in written["models"]] == [model_a.key, model_b.key]
    assert [entry["color"] for entry in written["models"]] == qualitative.Plotly[:2]
    assert "visible" not in written["models"][0]  # top-1 visible by default
    assert written["models"][1]["visible"] is False


def test_write_site_payload_subset_run_merges(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Subset runs splice fresh entries into the committed payload by id (replacing
    stale ones, appending new models, pruning retired ones), take shared fields fresh
    and sort entries into active-roster order by default.
    """
    model_a, model_b = list(Model.active())[:2]
    committed = {
        "shared": "old",
        "models": [
            {"key": model_b.key, "y": [2]},  # to be replaced by the fresh entry
            {"key": "retired-model", "y": [0]},  # to be pruned
            {"key": model_a.key, "y": [1]},  # to be kept as committed
        ],
    }
    figs.write_json_gz(f"{site_fig_dir}/demo.json.gz", committed)
    monkeypatch.setattr(cli_args, "models", [model_b])  # single-model run
    fresh = {"shared": "new", "models": [{"key": model_b.key, "y": [3]}]}
    figs.write_site_payload("demo", fresh)
    merged = load_payload(f"{site_fig_dir}/demo.json.gz")
    assert merged["shared"] == "new"  # shared fields are model-independent -> fresh
    assert merged["models"] == [  # sorted into active-roster order
        {"key": model_a.key, "y": [1]},
        {"key": model_b.key, "y": [3]},
    ]


def test_write_site_payload_subset_run_requires_committed_payload(
    site_fig_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Subset runs can't create a payload from scratch (would ship a partial roster)."""
    assert site_fig_dir.is_dir()
    monkeypatch.setattr(cli_args, "models", [next(iter(Model.active()))])
    with pytest.raises(FileNotFoundError, match="merge into an existing payload"):
        figs.write_site_payload("missing", {"models": []})
