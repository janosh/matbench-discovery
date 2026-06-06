"""Export analysis data as compact gzipped JSON payloads for the Svelte site.

Payloads are *data-only*: series arrays plus data-derived stats (MAE, AUC, F1, ...).
All presentation (axes, ref lines, legends, colors of fixed series) lives inline in
the Svelte pages that import these files. site/src/figs/payloads.d.ts documents the
expected payload shapes. Payloads are committed at site/src/figs/<name>.json.gz and
imported by pages through the json_gz Vite plugin (site/vite.config.ts), so each
route's chunk only contains the data of the figures it renders.

Helpers:
- ``write_json_gz(path, payload)``: deterministic gzipped JSON writer
- ``histogram(values)``: bin raw values into ``{x, y, bar_width}`` (HistBins shape)
- ``lttb`` / ``downsample_points``: down-sample over-resolved line/point series
- ``trace_xy`` / ``trace_color`` / ``trace_visible``: pull arrays/styles out of
  plotly traces (for figs built by pymatviz/matbench_discovery.plots helpers)
- ``sunburst_tree`` / ``sankey_data``: convert plotly sunburst/sankey figures into
  the nested/flat structures matterviz Sunburst/Sankey render directly
"""

from __future__ import annotations

import base64
import gzip
import json
import math
import os
from typing import TYPE_CHECKING, Any, Final

import numpy as np

if TYPE_CHECKING:
    import plotly.graph_objects as go

COORD_DECIMALS: Final = 5
DEFAULT_HIST_BINS: Final = 100


def round_list(values: Any) -> list[Any]:
    """Convert an array-like to a JSON-safe list: round floats to COORD_DECIMALS,
    keep ints/strings, replace non-finite numbers with None.
    """
    return [
        (round(val, COORD_DECIMALS) if math.isfinite(val) else None)
        if isinstance(val, float)
        else val
        for val in np.asarray(values).tolist()
    ]


def decode_array(value: Any) -> np.ndarray | None:
    """Decode a plotly value (plain list/array or base64 typed array) to numpy."""
    if value is None:
        return None
    if isinstance(value, dict) and "bdata" in value and "dtype" in value:
        raw = base64.b64decode(value["bdata"])
        dtype = np.dtype(value["dtype"]).newbyteorder("<")
        arr = np.frombuffer(raw, dtype=dtype)
        shape = value.get("shape")
        if shape is not None:
            dims = tuple(int(dim) for dim in str(shape).replace(" ", "").split(","))
            arr = arr.reshape(dims)
        return arr
    return np.asarray(value)


# === down-sampling ===
def lttb(x: np.ndarray, y: np.ndarray, n_out: int) -> tuple[np.ndarray, np.ndarray]:
    """Largest-Triangle-Three-Buckets down-sampling preserving visual line shape."""
    n = len(x)
    if n_out >= n or n_out < 3:
        return x, y
    sampled_idx = [0]
    bucket_size = (n - 2) / (n_out - 2)
    pos_a = 0
    for bucket in range(n_out - 2):
        lo = math.floor((bucket + 1) * bucket_size) + 1
        hi = min(math.floor((bucket + 2) * bucket_size) + 1, n)
        if hi > lo:
            avg_x, avg_y = float(np.mean(x[lo:hi])), float(np.mean(y[lo:hi]))
            candidates = range(lo, hi)
        else:
            avg_x, avg_y = float(x[-1]), float(y[-1])
            candidates = range(lo, lo + 1)
        ax, ay = float(x[pos_a]), float(y[pos_a])
        best_area, best_idx = -1.0, lo
        for cand in candidates:
            if cand >= n:
                break
            area = abs(
                (ax - avg_x) * (float(y[cand]) - ay)
                - (ax - float(x[cand])) * (avg_y - ay)
            )
            if area > best_area:
                best_area, best_idx = area, cand
        sampled_idx.append(best_idx)
        pos_a = best_idx
    sampled_idx.append(n - 1)
    idx = np.asarray(sorted(set(sampled_idx)))
    return x[idx], y[idx]


def downsample_points(
    x: np.ndarray, y: np.ndarray, n_out: int, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Randomly (seeded) down-sample an unordered point cloud to at most ``n_out``."""
    n = len(x)
    if n <= n_out:
        return x, y
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=n_out, replace=False))
    return x[idx], y[idx]


# === data builders ===
def histogram(
    values: Any,
    *,
    bins: int = DEFAULT_HIST_BINS,
    value_range: tuple[float, float] | None = None,
) -> dict[str, Any]:
    """Bin raw values into the HistBins payload shape: pre-computed bin centers,
    counts and bar width (drops NaN/inf values before binning).
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    counts, edges = np.histogram(arr, bins=bins, range=value_range)
    return {
        "x": round_list((edges[:-1] + edges[1:]) / 2),
        "y": [int(count) for count in counts],
        "bar_width": round(float(edges[1] - edges[0]), 6),
    }


# === plotly trace extraction (for pymatviz-built figures) ===
def trace_xy(trace: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return a plotly trace's x/y as numpy arrays (decoding typed arrays)."""
    x, y = decode_array(trace.x), decode_array(trace.y)
    if x is None or y is None:
        raise ValueError(f"trace {getattr(trace, 'name', '')!r} has no x/y data")
    return x, y


def trace_color(trace: Any) -> str | None:
    """Best-effort single color for a trace (line color, else marker color)."""
    line_color = getattr(getattr(trace, "line", None), "color", None)
    if isinstance(line_color, str):
        return line_color
    marker_color = getattr(getattr(trace, "marker", None), "color", None)
    return marker_color if isinstance(marker_color, str) else None


def trace_visible(trace: Any) -> bool:
    """Map plotly ``visible`` (True/False/'legendonly') to a boolean."""
    vis = getattr(trace, "visible", True)
    return vis is True or vis is None


def _get_trace(fig: go.Figure | dict[str, Any], trace_type: str) -> dict[str, Any]:
    """Find the first trace of ``trace_type`` in a plotly figure (or fig dict)."""
    fig_dict = fig if isinstance(fig, dict) else fig.to_plotly_json()
    trace = next(
        (tr for tr in fig_dict.get("data", []) if tr.get("type") == trace_type), None
    )
    if trace is None:
        raise ValueError(f"figure contains no {trace_type} trace")
    return trace


def sunburst_tree(fig: go.Figure | dict[str, Any]) -> list[dict[str, Any]]:
    """Convert a plotly sunburst figure (flat labels/parents/values) into the nested
    node list matterviz Sunburst renders (SunburstNode[]).
    """
    trace = _get_trace(fig, "sunburst")
    labels = [str(val) for val in round_list(decode_array(trace.get("labels")))]
    parents = [str(val) for val in round_list(decode_array(trace.get("parents")))]
    values = round_list(decode_array(trace.get("values")))
    raw_ids = trace.get("ids")
    ids = (
        [str(val) for val in round_list(decode_array(raw_ids))]
        if raw_ids is not None
        else labels
    )
    nodes: dict[str, dict[str, Any]] = {
        node_id: {
            "id": node_id,
            "label": labels[idx],
            "value": values[idx],
            "children": [],
        }
        for idx, node_id in enumerate(ids)
    }
    roots: list[dict[str, Any]] = []
    for idx, node_id in enumerate(ids):
        parent = parents[idx] if idx < len(parents) else ""
        siblings = nodes[parent]["children"] if parent in nodes else roots
        siblings.append(nodes[node_id])
    return roots


def sankey_data(fig: go.Figure | dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Convert a plotly sankey figure into matterviz Sankey {nodes, links}."""
    trace = _get_trace(fig, "sankey")
    node, link = trace.get("node") or {}, trace.get("link") or {}
    labels = round_list(decode_array(node.get("label")))
    nodes: list[dict[str, Any]] = [{"label": str(label)} for label in labels]
    colors = node.get("color")
    if isinstance(colors, (list, dict)):
        for box, color in zip(nodes, round_list(decode_array(colors)), strict=False):
            if isinstance(color, str):
                box["color"] = color
    sources = round_list(decode_array(link.get("source")))
    targets = round_list(decode_array(link.get("target")))
    values = round_list(decode_array(link.get("value")))
    links = [
        {"source": int(src), "target": int(tgt), "value": float(val)}
        for src, tgt, val in zip(sources, targets, values, strict=True)
    ]
    if not nodes or not links:
        raise ValueError("sankey trace has no nodes or links")
    return {"nodes": nodes, "links": links}


# === IO ===
def write_json_gz(path: str, data: dict[str, Any]) -> int:
    """Write deterministic gzipped JSON; return compressed byte size."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
    compressed = gzip.compress(payload, compresslevel=9, mtime=0)
    with open(path, "wb") as file:
        file.write(compressed)
    return len(compressed)
