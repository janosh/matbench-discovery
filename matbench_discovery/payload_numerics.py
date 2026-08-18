"""Deterministic numerical transforms used by generated site payloads."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any, Final

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping

    import numpy.typing as npt

COORD_DECIMALS: Final = 5


def round_list(values: npt.ArrayLike | None) -> list[Any]:
    """Convert values to a JSON-safe list with bounded floating precision."""
    if values is None:
        return []
    return [
        (round(value, COORD_DECIMALS) if math.isfinite(value) else None)
        if isinstance(value, float)
        else value
        for value in np.asarray(values).tolist()
    ]


def evenly_spaced_indices(n_values: int, n_out: int) -> np.ndarray:
    """Return at most ``n_out`` distinct indices spanning ``n_values`` entries."""
    if n_values < 1 or n_out < 1:
        raise ValueError(f"expected positive sizes, got {n_values=} and {n_out=}")
    return np.linspace(0, n_values - 1, min(n_values, n_out), dtype=int)


def lttb(
    x_values: np.ndarray, y_values: np.ndarray, n_out: int
) -> tuple[np.ndarray, np.ndarray]:
    """Down-sample a line with Largest-Triangle-Three-Buckets."""
    n_values = len(x_values)
    if n_out >= n_values or n_out < 3:
        return x_values, y_values
    sampled_indices = [0]
    bucket_size = (n_values - 2) / (n_out - 2)
    previous_idx = 0
    for bucket_idx in range(n_out - 2):
        average_start = math.floor((bucket_idx + 1) * bucket_size) + 1
        average_end = min(math.floor((bucket_idx + 2) * bucket_size) + 1, n_values)
        avg_x = float(np.mean(x_values[average_start:average_end]))
        avg_y = float(np.mean(y_values[average_start:average_end]))
        candidate_start = math.floor(bucket_idx * bucket_size) + 1
        candidate_end = math.floor((bucket_idx + 1) * bucket_size) + 1
        anchor_x = float(x_values[previous_idx])
        anchor_y = float(y_values[previous_idx])
        best_area, best_idx = -1.0, candidate_start
        for candidate_idx in range(candidate_start, candidate_end):
            area = abs(
                (anchor_x - avg_x) * (float(y_values[candidate_idx]) - anchor_y)
                - (anchor_x - float(x_values[candidate_idx])) * (avg_y - anchor_y)
            )
            if area > best_area:
                best_area, best_idx = area, candidate_idx
        sampled_indices.append(best_idx)
        previous_idx = best_idx
    sampled_indices.append(n_values - 1)
    indices = np.asarray(sampled_indices, dtype=int)
    return x_values[indices], y_values[indices]


def histogram(
    values: npt.ArrayLike,
    *,
    bins: int,
    value_range: tuple[float, float] | None = None,
) -> dict[str, Any]:
    """Bin finite values into centers, counts, and a bar width."""
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    counts, edges = np.histogram(array, bins=bins, range=value_range)
    return {
        "x": round_list((edges[:-1] + edges[1:]) / 2),
        "y": counts.tolist(),
        "bar_width": round(float(edges[1] - edges[0]), 6),
    }


def sankey_payload_from_flow(flow_data: Mapping[str, Any]) -> dict[str, Any]:
    """Canonicalize ``pymatviz.sankey_flow_data`` output for a site payload."""
    labels = [str(label) for label in flow_data["labels"]]
    source_indices = np.asarray(flow_data["source_indices"], dtype=int).tolist()
    target_indices = np.asarray(flow_data["target_indices"], dtype=int).tolist()
    values = round_list(flow_data["value"])
    if not labels or not source_indices or not target_indices:
        raise ValueError("sankey flow has no nodes or links")
    used_indices = sorted({*source_indices, *target_indices})
    remap = {old_idx: new_idx for new_idx, old_idx in enumerate(used_indices)}
    links = sorted(
        zip(
            (remap[source] for source in source_indices),
            (remap[target] for target in target_indices),
            values,
            strict=True,
        )
    )
    sources, targets, values = map(list, zip(*links, strict=True))
    return {
        "labels": [labels[idx] for idx in used_indices],
        "source": sources,
        "target": targets,
        "value": values,
    }
