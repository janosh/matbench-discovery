"""Export analysis data as compact gzipped JSON payloads for the Svelte site.

Payloads are *data-only*: series arrays plus data-derived stats (MAE, AUC, F1, ...).
All presentation (axes, ref lines, legends, per-model colors, render order, default
visibility) lives inline in the Svelte pages that import these files.
site/src/figs/payloads.d.ts documents the expected payload shapes.

Static payloads are committed as a single gzipped site/src/figs/<name>.json.gz.
Multi-model payloads are committed as line-delimited site/src/figs/<name>.jsonl: one
JSON object per line (a lone ``{"_base": {...}}`` line for shared fields + one line per
model). Two submissions that each add a model insert different lines, which git merges
cleanly rather than colliding on one un-mergeable gzipped blob. The figure_payload
Vite plugin (site/vite.config.ts) loads both; .jsonl reassembles into the aggregate.

Helpers:
- ``write_json_gz(path, payload)``: deterministic gzipped JSON writer
- ``write_site_payload(name, payload)``: write a multi-model payload as line-delimited
  JSONL (position-independent data; presentation applied client-side)
- ``read_jsonl_payload(path)``: reassemble a .jsonl payload into ``{**shared, models}``
- ``histogram(values)``: bin raw values into ``{x, y, bar_width}`` (HistBins shape)
- ``lttb``: down-sample over-resolved line series
- ``trace_xy`` / ``trace_color`` / ``trace_payload``: pull arrays/styles out of
  plotly traces (for figs built by pymatviz/plots helpers)
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
    import numpy.typing as npt
    import plotly.graph_objects as go
    from plotly.basedatatypes import BaseTraceType

COORD_DECIMALS: Final = 5
DEFAULT_HIST_BINS: Final = 100


def round_list(values: npt.ArrayLike | None) -> list[Any]:
    """Convert an array-like to a JSON-safe list: round floats to COORD_DECIMALS,
    keep ints/strings, replace non-finite numbers with None. ``None`` -> ``[]`` (so a
    missing trace field surfaces as an empty list, not a TypeError on iteration).

    Don't replace calls with plain ``.tolist()``: rounding cuts gzipped payload sizes
    to ~1/3 (full float64 repr keeps 17 significant digits) and NaN -> None is required
    since write_json_gz uses ``json.dumps(allow_nan=False)`` which raises on NaN.
    """
    if values is None:
        return []
    return [
        (round(val, COORD_DECIMALS) if math.isfinite(val) else None)
        if isinstance(val, float)
        else val
        for val in np.asarray(values).tolist()
    ]


def str_list(values: npt.ArrayLike | dict[str, str] | None) -> list[str]:
    """Convert a (possibly base64-encoded plotly) array of labels to a list of str.
    ``None`` -> ``[]``.
    """
    arr = decode_array(values)
    return [] if arr is None else [str(val) for val in arr.tolist()]


def decode_array(value: npt.ArrayLike | dict[str, str] | None) -> np.ndarray | None:
    """Decode a plotly value (plain list/array or base64 typed array) to numpy."""
    if value is None:
        return None
    if isinstance(value, dict) and "bdata" in value and "dtype" in value:
        # plotly base64 typed-array payload; coerce to a plain str dict so indexing is
        # typed (isinstance narrowing alone yields dict[Unknown] via npt.ArrayLike)
        meta = {str(key): str(val) for key, val in value.items()}
        raw = base64.b64decode(meta["bdata"])
        dtype = np.dtype(meta["dtype"]).newbyteorder("<")
        arr = np.frombuffer(raw, dtype=dtype)
        shape = meta.get("shape")
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
    values: npt.ArrayLike,
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
        "y": counts.tolist(),  # np.histogram counts are already int64
        "bar_width": round(float(edges[1] - edges[0]), 6),
    }


# === plotly trace extraction (for pymatviz-built figures) ===
def trace_xy(trace: BaseTraceType) -> tuple[np.ndarray, np.ndarray]:
    """Return a plotly trace's x/y as numpy arrays (decoding typed arrays)."""
    x = decode_array(getattr(trace, "x", None))
    y = decode_array(getattr(trace, "y", None))
    if x is None or y is None:
        raise ValueError(f"trace {getattr(trace, 'name', '')!r} has no x/y data")
    return x, y


def trace_color(trace: BaseTraceType) -> str | None:
    """Best-effort single color for a trace (line color, else marker color)."""
    line_color = getattr(getattr(trace, "line", None), "color", None)
    if isinstance(line_color, str):
        return line_color
    marker_color = getattr(getattr(trace, "marker", None), "color", None)
    return marker_color if isinstance(marker_color, str) else None


def trace_payload(trace: BaseTraceType, *, x: bool = True) -> dict[str, Any]:
    """Standard payload entry for a plotly trace: label, color (if any), x/y arrays.

    Pass ``x=False`` for payloads whose models share a single top-level x array.
    """
    x_arr, y_arr = trace_xy(trace)
    entry: dict[str, Any] = {"label": str(getattr(trace, "name", ""))}
    if color := trace_color(trace):
        entry["color"] = color
    if x:
        entry["x"] = round_list(x_arr)
    return entry | {"y": round_list(y_arr)}


def _get_trace(fig: go.Figure | dict[str, Any], trace_type: str) -> dict[str, Any]:
    """Find the first trace of ``trace_type`` in a plotly figure (or fig dict)."""
    fig_dict: dict[str, Any] = fig if isinstance(fig, dict) else fig.to_plotly_json()  # ty: ignore[invalid-assignment]
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
        "labels": str_list(trace.get("labels")),
        "parents": str_list(trace.get("parents")),
        "values": round_list(decode_array(trace.get("values"))),
    }
    if (raw_ids := trace.get("ids")) is not None:
        out["ids"] = str_list(raw_ids)
    return out


def sankey_data(fig: go.Figure | dict[str, Any]) -> dict[str, Any]:
    """Extract a plotly sankey's node labels + link source/target/value arrays.

    matterviz's sankey_from_links builds the ``{nodes, links}`` SankeyData from these
    flat arrays, so shipping them keeps payloads smaller + avoids duplicating its logic.
    """
    trace = _get_trace(fig, "sankey")
    node, link = trace.get("node") or {}, trace.get("link") or {}
    labels = str_list(node.get("label"))
    sources = decode_array(link.get("source"))
    targets = decode_array(link.get("target"))
    if not labels or sources is None or targets is None or len(sources) == 0:
        raise ValueError("sankey trace has no nodes or links")
    # link indices are ints; round_list would map a non-finite to None -> int(None)
    src_idx = np.asarray(sources, dtype=int).tolist()
    tgt_idx = np.asarray(targets, dtype=int).tolist()
    # drop nodes no link references (plotly keeps many unused spacegroup nodes whose
    # crammed labels overlap) and reindex the links onto the kept nodes
    used = sorted({*src_idx, *tgt_idx})
    remap = {old: new for new, old in enumerate(used)}
    # canonicalize link order: upstream pandas value_counts can permute equal-count
    # links between runs, which would make payload bytes depend on run composition
    links = sorted(
        zip(
            (remap[src] for src in src_idx),
            (remap[tgt] for tgt in tgt_idx),
            round_list(decode_array(link.get("value"))),
            strict=True,
        )
    )
    sources, targets, values = map(list, zip(*links, strict=True))
    return {
        "labels": [labels[idx] for idx in used],
        "source": sources,
        "target": targets,
        "value": values,
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


def read_jsonl_payload(path: str) -> dict[str, Any]:
    """Read a JSONL figure payload (from ``write_site_payload``) into the aggregate
    ``{**shared, models: [...]}`` dict. One JSON object per line; the lone
    ``{"_base": {...}}`` line carries shared fields, every other line is a model entry.
    Mirrors the site's jsonl Vite loader (site/vite.config.ts).
    """
    base: dict[str, Any] = {}
    models: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as file:
        for line in file:
            if not (line := line.strip()):
                continue
            entry = json.loads(line)
            if set(entry) == {"_base"}:
                base = entry["_base"]
            else:
                models.append(entry)
    return {**base, "models": models}


def write_jsonl_payload(
    path: str, payload: dict[str, Any], *, id_field: str = "key", full_run: bool
) -> int:
    """Write a multi-model payload as line-delimited JSONL at ``path`` (shared writer
    behind ``write_site_payload``; also used for the per-element-errors payload).

    One JSON object per line - a lone ``{"_base": {...}}`` shared-fields line plus one
    line per model, sorted by ``id_field`` and stripped of presentation (applied
    client-side). ``full_run`` rewrites the whole roster; subset --models runs splice
    fresh entries into the committed file by ``id_field``, keeping the committed _base.
    """

    def model_id(model: dict[str, Any]) -> str:
        return str(model.get(id_field) or model["label"])

    # strip presentation fields so lines stay position-independent (set client-side)
    models = [
        {key: val for key, val in model.items() if key not in ("color", "visible")}
        for model in payload["models"]
    ]
    shared_from = payload  # full runs take shared fields from the fresh payload
    if not full_run:  # splice fresh entries into the committed file by id
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"{path} not found: subset runs (--models) splice into an existing "
                "payload. Run without --models to regenerate from scratch."
            )
        committed = read_jsonl_payload(path)
        fresh = {model_id(model): model for model in models}
        models = [fresh.pop(model_id(old), old) for old in committed["models"]]
        models += list(fresh.values())
        # shared fields are model-independent; keep the committed _base so a subset run
        # only rewrites model lines (no churn from fresh dict order / float noise)
        shared_from = committed
    shared = {key: val for key, val in shared_from.items() if key != "models"}

    models.sort(key=model_id)
    # a lone _base line (when shared fields exist) followed by one line per model
    records = [{"_base": shared}, *models] if shared else models
    body = "".join(
        json.dumps(record, allow_nan=False, separators=(",", ":")) + "\n"
        for record in records
    )
    if dir_name := os.path.dirname(path):
        os.makedirs(dir_name, exist_ok=True)
    with open(path, "w", encoding="utf-8") as file:
        file.write(body)
    n_bytes = len(body.encode())
    print(f"Wrote {os.path.basename(path)} ({n_bytes:,} bytes, {len(models)} models)")
    return n_bytes


def write_site_payload(
    name: str, payload: dict[str, Any], *, id_field: str = "key"
) -> int:
    """Write a multi-model figure payload as one JSONL file, site/src/figs/<name>.jsonl
    (thin wrapper over ``write_jsonl_payload``; full vs subset run from CLI --models).
    """
    from matbench_discovery import SITE_FIG_DATA
    from matbench_discovery.cli import is_full_model_run

    path = f"{SITE_FIG_DATA}/{name}.jsonl"
    return write_jsonl_payload(
        path, payload, id_field=id_field, full_run=is_full_model_run()
    )
