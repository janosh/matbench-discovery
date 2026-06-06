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
- ``lttb``: down-sample over-resolved line series
- ``trace_xy`` / ``trace_color`` / ``trace_visible`` / ``trace_payload``: pull
  arrays/styles out of plotly traces (for figs built by pymatviz/plots helpers)
- ``sunburst_data`` / ``sankey_data``: pull the flat arrays out of plotly sunburst/
  sankey figures (matterviz builds the nested structures from these client-side)
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


def trace_payload(trace: Any, *, x: bool = True) -> dict[str, Any]:
    """Standard payload entry for a plotly trace: label, color (if any), x/y arrays.

    Pass ``x=False`` for payloads whose models share a single top-level x array.
    """
    x_arr, y_arr = trace_xy(trace)
    entry: dict[str, Any] = {"label": str(trace.name)}
    if color := trace_color(trace):
        entry["color"] = color
    if x:
        entry["x"] = round_list(x_arr)
    return entry | {"y": round_list(y_arr)}


def _get_trace(fig: go.Figure | dict[str, Any], trace_type: str) -> dict[str, Any]:
    """Find the first trace of ``trace_type`` in a plotly figure (or fig dict)."""
    fig_dict = fig if isinstance(fig, dict) else fig.to_plotly_json()
    trace = next(
        (tr for tr in fig_dict.get("data", []) if tr.get("type") == trace_type), None
    )
    if trace is None:
        raise ValueError(f"figure contains no {trace_type} trace")
    return trace


def sunburst_data(fig: go.Figure | dict[str, Any]) -> dict[str, Any]:
    """Extract a plotly sunburst's flat labels/parents/values(/ids) arrays.

    matterviz's ``sunburst_from_labels_parents`` builds the nested SunburstNode tree
    from these client-side (and handles duplicate ids), so shipping the flat arrays
    keeps payloads ~20% smaller than a pre-nested tree and avoids duplicating its logic.
    """
    trace = _get_trace(fig, "sunburst")
    out: dict[str, Any] = {
        "labels": [str(val) for val in round_list(decode_array(trace.get("labels")))],
        "parents": [str(val) for val in round_list(decode_array(trace.get("parents")))],
        "values": round_list(decode_array(trace.get("values"))),
    }
    if (raw_ids := trace.get("ids")) is not None:
        out["ids"] = [str(val) for val in round_list(decode_array(raw_ids))]
    return out


def sankey_data(fig: go.Figure | dict[str, Any]) -> dict[str, Any]:
    """Extract a plotly sankey's node labels + link source/target/value arrays.

    matterviz's sankey_from_links builds the {nodes, links} SankeyData from these flat
    arrays, so shipping them keeps payloads smaller and avoids duplicating its logic.
    """
    trace = _get_trace(fig, "sankey")
    node, link = trace.get("node") or {}, trace.get("link") or {}
    labels = [str(label) for label in round_list(decode_array(node.get("label")))]
    sources = round_list(decode_array(link.get("source")))
    targets = round_list(decode_array(link.get("target")))
    values = round_list(decode_array(link.get("value")))
    if not labels or not sources:
        raise ValueError("sankey trace has no nodes or links")
    return {
        "labels": labels,
        "source": [int(src) for src in sources],
        "target": [int(tgt) for tgt in targets],
        "value": [float(val) for val in values],
    }


# === IO ===
def write_json_gz(path: str, data: dict[str, Any]) -> int:
    """Write deterministic gzipped JSON; return compressed byte size."""
    if dir_name := os.path.dirname(path):  # empty for a bare filename -> skip makedirs
        os.makedirs(dir_name, exist_ok=True)
    payload = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
    compressed = gzip.compress(payload, compresslevel=9, mtime=0)
    with open(path, "wb") as file:
        file.write(compressed)
    return len(compressed)
