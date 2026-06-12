"""Export analysis data as compact gzipped JSON payloads for the Svelte site.

Payloads are *data-only*: series arrays plus data-derived stats (MAE, AUC, F1, ...).
All presentation (axes, ref lines, legends, colors of fixed series) lives inline in
the Svelte pages that import these files. site/src/figs/payloads.d.ts documents the
expected payload shapes. Payloads are committed at site/src/figs/<name>.json.gz and
imported by pages through the json_gz Vite plugin (site/vite.config.ts), so each
route's chunk only contains the data of the figures it renders.

Helpers:
- ``write_json_gz(path, payload)``: deterministic gzipped JSON writer
- ``write_site_payload(name, payload)``: write a multi-model payload; subset runs
  (--models) merge into the committed payload instead of clobbering it
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
    from collections.abc import Callable

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


def write_site_payload(
    name: str,
    payload: dict[str, Any],
    *,
    id_field: str = "key",
    sort_key: Callable[[dict[str, Any]], Any] | None = None,
    assign_colors: bool = False,
    visible_top_n: int | None = None,
) -> int:
    """Write a multi-model figure payload to site/src/figs/<name>.json.gz.

    Full runs (--models covers every active model) overwrite the payload wholesale.
    Subset runs (e.g. single-model ingestion) instead splice the freshly computed
    model entries into the committed payload by ``id_field``, so a contributor can
    refresh their own model without every other model's prediction files
    (https://github.com/janosh/matbench-discovery/issues/342). Entries of models that
    left the active roster (superseded) are pruned to keep the roster guards in
    tests/test_fig_payloads.py satisfiable; non-'models' top-level fields are
    model-independent reference data and taken fresh.

    Presentation fields are (re)assigned deterministically so subset merges stay
    byte-identical to full regens: ``sort_key`` orders entries (default:
    active-roster order, the order full runs iterate models in), ``assign_colors``
    cycles the plotly palette in that order, ``visible_top_n`` hides all but the
    first n entries by default.
    """
    from matbench_discovery import SITE_FIG_DATA
    from matbench_discovery.cli import is_full_model_run
    from matbench_discovery.enums import Model

    path = f"{SITE_FIG_DATA}/{name}.json.gz"
    models: list[dict[str, Any]] = list(payload["models"])
    if not is_full_model_run():
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"{path} not found: subset runs (--models) can only merge into an "
                "existing payload. Run without --models to regenerate from scratch."
            )
        with gzip.open(path) as file:
            committed = json.load(file)
        # fresh entries replace committed ones, new models append at the end
        fresh_by_id = {entry[id_field]: entry for entry in models}
        old = committed["models"]
        models = [fresh_by_id.pop(entry[id_field], entry) for entry in old]
        models += list(fresh_by_id.values())

    roster = {getattr(model, id_field): idx for idx, model in enumerate(Model.active())}
    models = [entry for entry in models if entry[id_field] in roster]
    models.sort(key=sort_key or (lambda entry: roster[entry[id_field]]))
    if assign_colors:
        from plotly.express.colors import qualitative

        for idx, entry in enumerate(models):
            entry["color"] = qualitative.Plotly[idx % len(qualitative.Plotly)]
    if visible_top_n is not None:
        for idx, entry in enumerate(models):
            if idx < visible_top_n:  # visible defaults to true when absent
                entry.pop("visible", None)
            else:
                entry["visible"] = False

    n_bytes = write_json_gz(path, payload | {"models": models})
    print(f"Wrote {name}.json.gz ({n_bytes:,} bytes, {len(models)} models)")
    return n_bytes
