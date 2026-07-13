"""Export analysis data as compact gzipped JSON payloads for the Svelte site.

Payloads are *data-only*: series arrays plus data-derived stats (MAE, AUC, F1, ...).
All presentation (axes, ref lines, legends, per-model colors, render order, default
visibility) lives inline in the Svelte pages that import these files.
site/src/figs/payloads.d.ts documents the expected payload shapes.

Static payloads are committed as a single gzipped site/src/figs/<name>.json.gz.
Multi-model payloads are committed as line-delimited site/src/figs/<name>.jsonl: one
JSON object per line (a lone ``{"_base": {...}}`` line for shared fields + one line per
model). Two submissions that each add a model insert different lines, which git merges
cleanly rather than colliding on one un-mergeable gzipped blob. The json_payload
Vite plugin (site/vite.config.ts) loads both; .jsonl reassembles into the aggregate.

Helpers:
- ``write_json_gz(path, payload)``: deterministic gzipped JSON writer
- ``write_site_payload(name, payload)``: write a multi-model payload as line-delimited
  JSONL (position-independent data; presentation applied client-side)
- ``read_jsonl_payload(path)``: reassemble a .jsonl payload into ``{**shared, models}``
- ``histogram(values)``: bin raw values into ``{x, y, bar_width}`` (HistBins shape)
- ``lttb``: down-sample over-resolved line series
- ``sankey_payload_from_flow``: canonicalize ``pymatviz.sankey_flow_data`` output
"""

from __future__ import annotations

import contextlib
import gzip
import json
import math
import os
import zlib
from typing import TYPE_CHECKING, Any, Final

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping

    import numpy.typing as npt

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


def sankey_payload_from_flow(flow_data: Mapping[str, Any]) -> dict[str, Any]:
    """Canonicalize ``pymatviz.sankey_flow_data`` output for a site payload."""
    labels = [str(label) for label in flow_data["labels"]]
    src_idx = np.asarray(flow_data["source_indices"], dtype=int).tolist()
    tgt_idx = np.asarray(flow_data["target_indices"], dtype=int).tolist()
    values = round_list(flow_data["value"])
    if not labels or not src_idx or not tgt_idx:
        raise ValueError("sankey flow has no nodes or links")
    used = sorted({*src_idx, *tgt_idx})
    remap = {old: new for new, old in enumerate(used)}
    links = sorted(
        zip(
            (remap[src] for src in src_idx),
            (remap[tgt] for tgt in tgt_idx),
            values,
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
    """Write deterministic gzipped JSON; return compressed byte size.

    Skips the write when the existing file already decompresses to the same JSON:
    different zlib builds encode identical input to different (valid) gzip streams,
    so unconditional rewrites would churn committed payloads across CI runners.
    """
    payload = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
    # a missing/corrupt existing file falls through to the (re)write below. Per the
    # gzip docs, invalid files raise OSError (BadGzipFile), EOFError (truncation)
    # or zlib.error (corrupt deflate stream)
    with contextlib.suppress(OSError, EOFError, zlib.error), gzip.open(path) as file:
        if file.read() == payload:
            return os.path.getsize(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)  # "." for bare filenames
    compressed = gzip.compress(payload, compresslevel=9, mtime=0)
    with open(path, "wb") as file:
        return file.write(compressed)


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
    fresh entries into the committed file by ``id_field``, keeping the committed _base
    and dropping entries of models that are no longer active (e.g. superseded by the
    spliced-in model - a splice alone could never prune them).
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
        # prune committed entries of known-but-inactive models (matched by key or
        # label so both id spaces are covered); unknown ids such as reference lines
        # ('Test set standard deviation') are preserved as-is
        from matbench_discovery.enums import Model

        inactive_ids = {
            id_attr
            for model in Model
            if not model.is_active
            for id_attr in (model.key, model.label)
        }
        models = [
            fresh.pop(model_id(old), old)
            for old in committed["models"]
            if model_id(old) not in inactive_ids
        ]
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
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
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
