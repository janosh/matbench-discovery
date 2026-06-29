"""Metrics for analyzing diatomic potential energy curves.

Big props to Tamas Stenczel and Yuan Chiang who spearheaded this type of PES
smoothness analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP
paper https://arxiv.org/abs/2401.00096 (see fig. 56) and MLIP Arena
https://huggingface.co/spaces/atomind/mlip-arena, respectively.
"""

from dataclasses import dataclass, field
from typing import Any, Self

import numpy as np
from ase.data import atomic_numbers, covalent_radii, vdw_alvarez
from numpy.typing import ArrayLike

from matbench_discovery.data import update_yaml_file
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.diatomics.energy import (
    calc_curve_diff_auc,
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae,
    calc_tortuosity,
)
from matbench_discovery.metrics.diatomics.force import (
    calc_conservation_deviation,
    calc_force_flips,
    calc_force_jump,
    calc_force_mae,
    calc_force_total_variation,
)


class DiatomicCurve:
    """Energies and forces for a single diatomic molecule at multiple distances."""

    distances: np.ndarray  # shape (n_distances,)
    energies: np.ndarray  # shape (n_distances,)
    forces: np.ndarray  # shape (n_distances, n_atoms, 3)

    def __init__(
        self, distances: ArrayLike, energies: ArrayLike, forces: ArrayLike
    ) -> None:
        """Convert inputs to numpy arrays and reshape forces if needed."""
        self.distances = np.asarray(distances)
        self.energies = np.asarray(energies)
        self.forces = np.asarray(forces)

        n_distances = len(self.distances)

        # Handle forces stored as (1, n_distances*n_atoms, 3)
        # instead of (n_distances, n_atoms, 3)
        if self.forces.shape == (1, 2 * n_distances, 3):
            self.forces = self.forces.reshape(n_distances, 2, 3)


@dataclass
class DiatomicCurves:
    """Container for diatomic potential energy curves and forces of multiple
    element pairs.

    Attributes:
        distances (np.ndarray): Interatomic distances in Å.
        homo_nuclear (dict[str, DiatomicCurve]): Map of element pairs
            (e.g. "H-H") to their DiatomicCurve (energies and forces).
        hetero_nuclear (dict[str, DiatomicCurve] | None): Optional map of element pairs
            (e.g. "H-He") to their DiatomicCurve.
    """

    distances: np.ndarray  # shape (n_distances,)
    homo_nuclear: dict[str, DiatomicCurve]
    hetero_nuclear: dict[str, DiatomicCurve] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create DiatomicCurves from a dictionary loaded from JSON."""
        dists = data["distances"] = np.asarray(data["distances"])
        for key in {"homo-nuclear", "hetero-nuclear"} & set(data):
            data[key.replace("-", "_")] = {
                formula: DiatomicCurve(
                    energies=dct["energies"], forces=dct["forces"], distances=dists
                )
                for formula, dct in data.pop(key).items()
                if len(dct["energies"]) > 0
            }
        return cls(**data)


def _eval_window(elem_symbol: str, seps_max: float) -> tuple[float, float]:
    """Element-specific [r_min, r_max] metric-evaluation window, matching MLIP Arena.

    MLIP Arena samples (and thus evaluates) homonuclear curves only over
    ``[0.9 * covalent_radius, 3.1 * vdW_radius]`` (capped at ``seps_max``), deliberately
    excluding the deep-overlap region (below ~the covalent radius) that no model saw in
    training and where energies/forces extrapolate to unphysical magnitudes. Without it
    every magnitude metric (energy/force jumps, grad-norm, total variation) is dominated
    by the steep repulsive wall at sub-Angstrom separations.

    Args:
        elem_symbol (str): Homonuclear pair label, e.g. "H-H".
        seps_max (float): Largest available separation, used to cap r_max.

    Returns:
        tuple[float, float]: (r_min, r_max) in Å.
    """
    atomic_num = atomic_numbers[elem_symbol.split("-", maxsplit=1)[0]]
    r_min = 0.9 * covalent_radii[atomic_num]
    r_vdw = (
        vdw_alvarez.vdw_radii[atomic_num]
        if atomic_num < len(vdw_alvarez.vdw_radii)
        else np.nan
    )
    r_max = 3.1 * r_vdw if np.isfinite(r_vdw) else seps_max
    return r_min, min(r_max, seps_max)


def calc_diatomic_metrics(
    ref_curves: DiatomicCurves | None,
    pred_curves: DiatomicCurves,
    metrics: dict[str, dict[str, Any]] | None = None,
    *,
    interpolate: bool | int = False,
) -> dict[str, dict[str, float]]:
    """Calculate diatomic curve metrics comparing predicted curves to reference curves.

    Args:
        ref_curves (DiatomicCurves | None): Reference energy curves for each element.
            If None, only metrics that don't require reference data will be calculated.
        pred_curves (DiatomicCurves): Predicted energy curves for each element.
        metrics (dict[str, dict[str, Any]] | None): Map of metric names to
            dictionaries of keyword arguments for each metric function. If None, uses
            all metrics with default parameters. To use a subset of metrics, provide
            a dictionary with those metric names as keys and their keyword arguments
            as values. Empty dictionaries will use default parameters.
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses default number of points for interpolation.
            If an integer, uses that many points for interpolation.
            This value is passed to each metric function that supports it.

    Returns:
        dict[str, dict[str, float]]: Map of element symbols to metric dicts with keys
            being the metric names and values being the metric values.
    """
    requested_metrics = set(metrics) if metrics is not None else set(MbdKey)
    if unknown_metrics := requested_metrics - set(MbdKey):
        raise ValueError(f"{unknown_metrics=}. Valid metrics=")

    results: dict[str, dict[str, float]] = {}
    metric_kwargs = (metrics or {}).copy()

    # Initialize empty kwargs for each metric if not provided
    for key in MbdKey:
        metric_kwargs.setdefault(key, {})

        # set user-specified interpolate on metrics that support it
        if key in (
            MbdKey.norm_auc,
            MbdKey.energy_mae,
            MbdKey.force_mae,
            MbdKey.conservation,
        ):
            metric_kwargs[key].setdefault("interpolate", interpolate)

    for elem_symbol, pred_data in pred_curves.homo_nuclear.items():
        # restrict to MLIP Arena's element-specific evaluation window first, excluding
        # the unphysical deep-overlap region that otherwise dominates every magnitude
        # metric (and where unstable models tend to produce non-finite values; checking
        # finiteness below the window means a curve only bad in deep overlap is kept)
        r_min, r_max = _eval_window(elem_symbol, float(np.max(pred_curves.distances)))
        pred_mask = (pred_curves.distances >= r_min) & (pred_curves.distances <= r_max)
        if int(pred_mask.sum()) < 5:  # too few points in window for stable metrics
            print(f"Skipping {elem_symbol} diatomic metrics: <5 points in eval window")
            continue
        seps_pred = pred_curves.distances[pred_mask]
        pred_energies = np.asarray(pred_data.energies)[pred_mask]
        pred_forces = np.asarray(pred_data.forces)[pred_mask]

        # skip curves still non-finite within the window (model instability at scored
        # geometries); they'd fail metric validation and abort the whole model otherwise
        if not (np.isfinite(pred_energies).all() and np.isfinite(pred_forces).all()):
            print(f"Skipping {elem_symbol} diatomic metrics: non-finite curve values")
            continue

        # (metric_key, function, positional args) for metrics needing only the
        # predicted curve
        metric_calls: list[tuple[MbdKey, Any, tuple[Any, ...]]] = [
            (MbdKey.tortuosity, calc_tortuosity, (seps_pred, pred_energies)),
            (
                MbdKey.energy_diff_flips,
                calc_energy_diff_flips,
                (seps_pred, pred_energies),
            ),
            (
                MbdKey.energy_grad_norm_max,
                calc_energy_grad_norm_max,
                (seps_pred, pred_energies),
            ),
            (MbdKey.energy_jump, calc_energy_jump, (seps_pred, pred_energies)),
            (
                MbdKey.conservation,
                calc_conservation_deviation,
                (seps_pred, pred_energies, pred_forces),
            ),
            (MbdKey.force_flips, calc_force_flips, (seps_pred, pred_forces)),
            (
                MbdKey.force_total_variation,
                calc_force_total_variation,
                (seps_pred, pred_forces),
            ),
            (MbdKey.force_jump, calc_force_jump, (seps_pred, pred_forces)),
        ]

        # prepend reference-requiring metrics when a matching reference curve exists
        if ref_curves and (ref_data := ref_curves.homo_nuclear.get(elem_symbol)):
            ref_mask = (ref_curves.distances >= r_min) & (ref_curves.distances <= r_max)
            seps_ref = ref_curves.distances[ref_mask]
            ref_energies = np.asarray(ref_data.energies)[ref_mask]
            ref_forces = np.asarray(ref_data.forces)[ref_mask]
            if not np.array_equal(seps_pred, seps_ref) and not interpolate:
                raise ValueError(
                    f"Reference and predicted distances must be same when "
                    f"{interpolate=}\n{seps_pred=}, {seps_ref=}"
                )
            metric_calls[:0] = [
                (
                    MbdKey.norm_auc,
                    calc_curve_diff_auc,
                    (seps_ref, ref_energies, seps_pred, pred_energies),
                ),
                (
                    MbdKey.energy_mae,
                    calc_energy_mae,
                    (seps_ref, ref_energies, seps_pred, pred_energies),
                ),
                (
                    MbdKey.force_mae,
                    calc_force_mae,
                    (seps_ref, ref_forces, seps_pred, pred_forces),
                ),
            ]

        results[elem_symbol] = {
            key: calc_fn(*args, **metric_kwargs[key])
            for key, calc_fn, args in metric_calls
            if key in requested_metrics
        }

    return results


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, dict[str, float]],
    pred_file_path: str | None = None,
    run_metadata: dict[str, str | float] | None = None,
) -> dict[str, str | float | None]:
    """Write diatomic metrics to model YAML file.

    Args:
        model (Model): Model to write metrics for.
        metrics (dict[str, dict[str, float]]): Map of element symbols to dicts of
            metric values.
        pred_file_path (str | None): If given, record this (repo-relative) path as the
            metrics.diatomics.pred_file. Otherwise an existing pred_file is preserved.
        run_metadata (dict[str, str | float] | None): Extra non-metric fields describing
            the prediction run (e.g. hardware, run_time_sec). Recorded ahead of the
            metric values; a recompute without run_metadata preserves existing values.

    Returns:
        dict[str, str | float | None]: The metrics.diatomics block written (file refs
            and run metadata first, then the mean of each metric across all elements).
    """
    from matbench_discovery import ROOT

    # mean of each metric over the elements that have a finite value: skips elements
    # whose windowed curve is degenerate (e.g. tortuosity is NaN for a flat curve), and
    # drops a metric entirely if no element has a finite value rather than writing an
    # invalid `.nan`. Union the keys since elements can differ (only some have a ref).
    metric_keys = list(dict.fromkeys(key for em in metrics.values() for key in em))
    mean_metrics: dict[str, str | float | None] = {}
    for metric in metric_keys:
        finite = [
            val
            for em in metrics.values()
            if (val := em.get(metric)) is not None and np.isfinite(val)
        ]
        if finite:
            mean_metrics[str(metric)] = float(f"{np.mean(finite):.4}")

    # carry over only recognized run metadata (it describes the source run, not the
    # computed metrics, so it stays valid on recalculation)
    existing = model.metrics.get("diatomics", {})
    existing = existing if isinstance(existing, dict) else {}
    run_metadata = run_metadata or {}
    pred_file_url = run_metadata.get("pred_file_url", existing.get("pred_file_url"))
    block: dict[str, str | float | None] = {}
    if pred_file_url is not None:
        if pred_file_path is not None:  # a passed path overrides any existing one
            block["pred_file"] = pred_file_path.removeprefix(f"{ROOT}/")
        elif "pred_file" in existing:
            block["pred_file"] = existing["pred_file"]
        block["pred_file_url"] = pred_file_url
    for key in ("hardware", "run_time_sec"):
        if key in existing:
            block[key] = existing[key]
    block |= run_metadata
    if "pred_file" in block and "pred_file_url" not in block:
        block.pop("pred_file")
    block |= mean_metrics

    # preserve_existing=False so a recompute fully replaces the block, dropping
    # deprecated metrics left in the YAML. This still runs when no finite metrics were
    # produced, preventing stale metric values from surviving.
    update_yaml_file(
        model.yaml_path, "metrics.diatomics", block, preserve_existing=False
    )
    print(f"Wrote {model.label} diatomic metrics to {model.yaml_path}")
    return block
