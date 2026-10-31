"""Reproducibility checks for registered task results.

Independent reruns of a task are scored with that task's own metric code and
compared to the scores recomputed from the section's registered ``pred_file`` (never
the YAML values). The verdict in ``metrics.<task>.reproducibility`` is ``reproduced``
iff every verdict metric changed by at most TOLERANCE (relative), ``drifted`` if
more, ``inconclusive`` if the registered file or a run lacks one. All diatomic
metrics decide the verdict. MD trajectory metrics (RDF/ADF/vDOS/pressure) are chaotic
under any perturbation of a 20 ps rollout, so only the frame-wise private-label
diagnostics (energy_rmse, force_rmse) decide and the rest are only listed.
"""

import gzip
import json
from collections.abc import Mapping, Sequence
from itertools import permutations
from typing import Any, Final, Literal

import pandas as pd

from matbench_discovery import today
from matbench_discovery.enums import Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.exclusions import drop_metric_exclusions

ReproTask = Literal["diatomics", "md"]
REPRO_CAUSES: Final = (
    "protocol_change",
    "environment_change",
    "nondeterministic",
    "unknown",
)
TOLERANCE: Final = 1e-3  # max relative change of a reproduced verdict metric
MD_FRAMEWISE_KEYS: Final = ("energy_rmse", "force_rmse")


def score_files(
    task: ReproTask, model: Model, paths: Sequence[str]
) -> list[dict[str, float]]:
    """Model-level scores of prediction files, computed like the published ones.

    Diatomics: PBE reference, interpolate=200, the model's curated exclusions
    dropped, finite means (scripts/evals/diatomic_metrics.py). MD: calc_md_metrics
    minus cost and n_systems over combined per-system CSVs, each covering exactly the
    systems of ``paths[0]`` (the registered file).
    """
    if task == "diatomics":
        ref_curves = diatomics.load_dft_reference_curves("PBE")
        scores = []
        for path in paths:
            with gzip.open(path, mode="rb") as file:
                pred = DiatomicCurves.from_dict(json.load(file))
            metrics = diatomics.calc_diatomic_metrics(ref_curves, pred, interpolate=200)
            metrics = drop_metric_exclusions(model.name, metrics)
            scores.append(diatomics.diatomic_metric_means(metrics))
        return scores
    frames = [
        pd.read_csv(path, index_col="system", float_precision="round_trip")
        for path in paths
    ]
    for path, df_md in zip(paths, frames, strict=True):
        if sorted(df_md.index) != (systems := sorted(frames[0].index)):
            raise ValueError(f"{path} covers {sorted(df_md.index)}, expected {systems}")
    return [
        {
            key: value
            for key, value in md_metrics.calc_md_metrics(df_md).items()
            if key in md_metrics.PER_SYSTEM_METRIC_COLS and isinstance(value, float)
        }
        for df_md in frames
    ]


def rel_change(value: float, reference: float) -> float:
    """``|value - reference| / |reference|``; a nonzero value vs a zero reference
    has no relative scale and counts as a 100% change (1.0).
    """
    if reference == 0:
        return float(value != 0)
    return abs(value - reference) / abs(reference)


def make_block(
    task: ReproTask,
    registered: Mapping[str, float],
    run_scores: Sequence[Mapping[str, float]],
    *,
    runs: Sequence[dict[str, Any]],
    repeat_scores: Sequence[Mapping[str, float]] = (),
    cause: str | None = None,
    note: str | None = None,
) -> dict[str, Any]:
    """Verdict block in YAML key order: per-metric max relative change over runs (3
    significant figures, zeros omitted) and, given >= 2 identical same-env repeats,
    the max change of the verdict metrics over every ordered repeat pair (noise
    floor). A cause is required iff drifted, a note unless reproduced.
    """
    rel_changes: dict[str, float] = {}
    for run in run_scores:
        for metric in registered.keys() & run.keys():
            change = rel_change(run[metric], registered[metric])
            rel_changes[metric] = max(rel_changes.get(metric, 0.0), change)
    all_scores = (registered, *run_scores)
    keys = MD_FRAMEWISE_KEYS if task == "md" else {k for s in all_scores for k in s}
    if any(key not in scores for key in keys for scores in all_scores):
        status = "inconclusive"
    else:
        drifted = any(rel_changes[key] > TOLERANCE for key in keys)
        status = "drifted" if drifted else "reproduced"
    needs = {"drifted": "a cause and a note", "inconclusive": "a note"}.get(status)
    if (cause is not None) != (status == "drifted") or (needs and not note):
        raise ValueError(
            f"{status} verdict needs {needs or 'no cause'}, got {cause=}, {note=} "
            f"(rel changes {rel_changes})"
        )
    block = {
        "checked": today,
        "status": status,
        "cause": cause,
        "note": note or None,
        "tolerance": TOLERANCE,
        "runs": list(runs),
        "rel_changes": {
            key: float(f"{val:.3g}") for key, val in rel_changes.items() if val
        },
    }
    if repeat_scores:
        repeat_change = max(
            rel_change(run_a[key], run_b[key])
            for run_a, run_b in permutations(repeat_scores, 2)
            for key in keys
        )
        block["repeat_rel_change"] = float(f"{repeat_change:.3g}")
    return {key: value for key, value in block.items() if value is not None}
