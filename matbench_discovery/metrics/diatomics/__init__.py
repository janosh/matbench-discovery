"""Metrics for analyzing diatomic potential energy curves.

Big props to Tamas Stenczel and Yuan Chiang who spearheaded this type of PES
smoothness analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP
paper https://arxiv.org/abs/2401.00096 (see fig. 56) and MLIP Arena
https://huggingface.co/spaces/atomind/mlip-arena, respectively.
"""

import gzip
import json
from dataclasses import dataclass, field
from typing import Any, Self

import numpy as np
from ase.data import atomic_numbers, covalent_radii, vdw_alvarez
from numpy.typing import ArrayLike

from matbench_discovery import repo_relative_path
from matbench_discovery.data import update_yaml_file
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.diatomics.energy import (
    calc_energy_diff_flips,
    calc_energy_jump,
    calc_pbe_bond_length_error,
    calc_pbe_energy_mae,
    calc_pbe_vib_freq_error,
    calc_pbe_wall_dist_mae,
    calc_pbe_well_depth_error,
    calc_tortuosity,
)
from matbench_discovery.metrics.diatomics.force import (
    calc_force_flips,
    calc_force_jump,
    calc_force_mae,
    calc_force_total_variation,
)

DiatomicsYamlValue = str | float | list[str] | None
DIATOMIC_METRIC_KEYS = frozenset(
    str(key)
    for key in (
        MbdKey.tortuosity,
        MbdKey.force_flips,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.force_total_variation,
        MbdKey.force_jump,
        MbdKey.pbe_wall_dist_mae,
        MbdKey.pbe_energy_mae,
        MbdKey.pbe_bond_length_error,
        MbdKey.pbe_well_depth_error,
        MbdKey.pbe_force_mae,
        MbdKey.pbe_vib_freq_error,
    )
)


def _homo_key(formula: str) -> str:
    """Return element key for homonuclear pair labels, else formula unchanged."""
    elem1, sep, elem2 = formula.partition("-")
    return elem1 if sep and elem1 == elem2 else formula


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
        distances = np.asarray(data["distances"])

        def make_curves(section: str) -> dict[str, DiatomicCurve]:
            """Convert one JSON curve section to DiatomicCurve objects."""
            curves = data.get(section, data.get(section.replace("-", "_"), {}))
            key_fn = _homo_key if section.startswith("homo") else str
            return {
                key_fn(formula): DiatomicCurve(
                    distances=distances,
                    energies=curve["energies"],
                    forces=curve.get("forces", []),
                )
                for formula, curve in curves.items()
                if len(curve["energies"]) > 0
            }

        return cls(
            distances=distances,
            homo_nuclear=make_curves("homo-nuclear"),
            hetero_nuclear=make_curves("hetero-nuclear"),
        )


def load_dft_reference_curves(
    functional: str = "PBE", ref_path: str | None = None
) -> DiatomicCurves:
    """Load bundled DFT reference curves as DiatomicCurves keyed by formula."""
    from matbench_discovery import ROOT

    dft_ref_path = ref_path or f"{ROOT}/site/src/lib/diatomics-dft.json.gz"
    open_fn = gzip.open if dft_ref_path.endswith(".gz") else open
    with open_fn(dft_ref_path, mode="rt", encoding="utf-8") as file:
        references = json.load(file)[functional]
    return DiatomicCurves(
        distances=np.array([]),
        homo_nuclear={
            _homo_key(formula): DiatomicCurve(
                distances=curve["distances"],
                energies=curve["energies"],
                forces=curve.get("forces", []),
            )
            for formula, curve in references.items()
        },
    )


def _eval_window(elem_symbol: str, seps_max: float) -> tuple[float, float]:
    """Element-specific [r_min, r_max] metric-evaluation window, matching MLIP Arena.

    MLIP Arena samples (and thus evaluates) homonuclear curves only over
    ``[0.9 * covalent_radius, 3.1 * vdW_radius]`` (capped at ``seps_max``), deliberately
    excluding the deep-overlap region (below ~the covalent radius) that no model saw in
    training and where energies/forces extrapolate to unphysical magnitudes. Without it
    every magnitude metric (energy/force jumps, total variation) is dominated
    by the steep repulsive wall at sub-Angstrom separations.

    Args:
        elem_symbol (str): Homonuclear pair label, e.g. "H-H".
        seps_max (float): Largest available separation, used to cap r_max.

    Returns:
        tuple[float, float]: (r_min, r_max) in Å.
    """
    atomic_num = atomic_numbers[elem_symbol.split("-", maxsplit=1)[0]]
    r_cov = covalent_radii[atomic_num] if atomic_num < len(covalent_radii) else np.nan
    r_min = 0.9 * r_cov if np.isfinite(r_cov) else 0.0
    vdw_radii = vdw_alvarez.vdw_radii
    r_vdw = vdw_radii[atomic_num] if atomic_num < len(vdw_radii) else np.nan
    r_max = min(3.1 * r_vdw, seps_max) if np.isfinite(r_vdw) else seps_max
    return r_min, r_max


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
    requested_metric_keys = (
        {str(key) for key in metrics} if metrics is not None else DIATOMIC_METRIC_KEYS
    )
    if unknown_metrics := requested_metric_keys - DIATOMIC_METRIC_KEYS:
        raise ValueError(f"{unknown_metrics=}. Valid metrics=")
    requested_metrics = {MbdKey(key) for key in requested_metric_keys}

    results: dict[str, dict[str, float]] = {}
    metric_kwargs = {
        MbdKey(key): kwargs.copy() for key, kwargs in (metrics or {}).items()
    }
    for key in requested_metrics & {
        MbdKey.pbe_energy_mae,
        MbdKey.pbe_force_mae,
    }:
        metric_kwargs.setdefault(key, {}).setdefault("interpolate", interpolate)

    for elem_symbol, pred_data in pred_curves.homo_nuclear.items():
        # restrict to MLIP Arena's element-specific evaluation window first, excluding
        # the unphysical deep-overlap region that otherwise dominates every magnitude
        # metric (and where unstable models tend to produce non-finite values; checking
        # finiteness below the window means a curve only bad in deep overlap is kept)
        pred_dists = pred_data.distances
        r_min, r_max = _eval_window(elem_symbol, float(np.max(pred_dists)))
        pred_mask = (pred_dists >= r_min) & (pred_dists <= r_max)
        if int(pred_mask.sum()) < 5:  # too few points in window for stable metrics
            print(f"Skipping {elem_symbol} diatomic metrics: <5 points in eval window")
            continue
        seps_pred = pred_dists[pred_mask]
        pred_energies = np.asarray(pred_data.energies)[pred_mask]
        pred_forces = np.asarray(pred_data.forces)[pred_mask]

        # skip curves still non-finite within the window (model instability at scored
        # geometries); they'd fail metric validation and abort the whole model otherwise
        if not (np.isfinite(pred_energies).all() and np.isfinite(pred_forces).all()):
            print(f"Skipping {elem_symbol} diatomic metrics: non-finite curve values")
            continue

        energy_args = (seps_pred, pred_energies)
        force_args = (seps_pred, pred_forces)
        # (metric_key, function, positional args) for metrics needing only the
        # predicted curve
        metric_calls: list[tuple[MbdKey, Any, tuple[Any, ...]]] = [
            (MbdKey.tortuosity, calc_tortuosity, energy_args),
            (
                MbdKey.energy_diff_flips,
                calc_energy_diff_flips,
                energy_args,
            ),
            (MbdKey.energy_jump, calc_energy_jump, energy_args),
            (MbdKey.force_flips, calc_force_flips, force_args),
            (
                MbdKey.force_total_variation,
                calc_force_total_variation,
                force_args,
            ),
            (MbdKey.force_jump, calc_force_jump, force_args),
        ]

        # prepend reference-requiring metrics when a matching reference curve exists
        if ref_curves and (ref_data := ref_curves.homo_nuclear.get(elem_symbol)):
            ref_dists = ref_data.distances
            ref_mask = (ref_dists >= r_min) & (ref_dists <= r_max)
            seps_ref = ref_dists[ref_mask]
            ref_energies = np.asarray(ref_data.energies)[ref_mask]
            if len(seps_ref) >= 2:
                energy_pair_args = (seps_ref, ref_energies, seps_pred, pred_energies)
                metric_calls[:0] = [
                    (metric_key, calc_fn, energy_pair_args)
                    for metric_key, calc_fn in (
                        (MbdKey.pbe_wall_dist_mae, calc_pbe_wall_dist_mae),
                        (MbdKey.pbe_energy_mae, calc_pbe_energy_mae),
                        (MbdKey.pbe_bond_length_error, calc_pbe_bond_length_error),
                        (MbdKey.pbe_well_depth_error, calc_pbe_well_depth_error),
                    )
                ] + [
                    (
                        MbdKey.pbe_vib_freq_error,
                        calc_pbe_vib_freq_error,
                        (elem_symbol, seps_ref, ref_energies, seps_pred, pred_energies),
                    ),
                ]
            ref_forces = np.asarray(ref_data.forces)
            if (
                ref_forces.size
                and len(ref_forces) == len(ref_dists)
                and len(seps_ref) >= 2
            ):
                ref_forces = ref_forces[ref_mask]
                force_interpolate = metric_kwargs.get(MbdKey.pbe_force_mae, {}).get(
                    "interpolate", False
                )
                same_grid = np.array_equal(seps_ref, seps_pred)
                has_overlap = max(seps_ref.min(), seps_pred.min()) <= min(
                    seps_ref.max(), seps_pred.max()
                )
                if same_grid or (force_interpolate and has_overlap):
                    force_pair_args = (seps_ref, ref_forces, seps_pred, pred_forces)
                    metric_calls.append(
                        (MbdKey.pbe_force_mae, calc_force_mae, force_pair_args)
                    )

        results[elem_symbol] = {
            key: calc_fn(*args, **metric_kwargs.get(key, {}))
            for key, calc_fn, args in metric_calls
            if key in requested_metrics
        }

    return results


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, dict[str, float]],
    pred_file_path: str | None = None,
    run_metadata: dict[str, str | float | list[str]] | None = None,
) -> dict[str, DiatomicsYamlValue]:
    """Write diatomic metrics to model YAML file.

    Args:
        model (Model): Model to write metrics for.
        metrics (dict[str, dict[str, float]]): Map of element symbols to dicts of
            metric values.
        pred_file_path (str | None): If given, record this path as
            metrics.diatomics.pred_file. Absolute paths must be inside the repo and are
            converted to repo-relative paths. Otherwise an existing pred_file is
            preserved.
        run_metadata (dict[str, str | float] | None): Extra non-metric fields describing
            the prediction run (e.g. hardware, run_time_sec). Recorded ahead of the
            metric values; a recompute without run_metadata preserves existing values.

    Returns:
        dict[str, DiatomicsYamlValue]: The metrics.diatomics block written (file refs
            and run metadata first, then metric means across all elements).
    """
    # mean of each metric over the elements that have a finite value: skips elements
    # whose windowed curve is degenerate (e.g. tortuosity is NaN for a flat curve), and
    # drops a metric entirely if no element has a finite value rather than writing an
    # invalid `.nan`. Union the keys since elements can differ (only some have a ref).
    mean_metrics: dict[str, DiatomicsYamlValue] = {}
    for metric in dict.fromkeys(
        key for elem_metrics in metrics.values() for key in elem_metrics
    ):
        finite = [
            val
            for elem_metrics in metrics.values()
            if (val := elem_metrics.get(metric)) is not None and np.isfinite(val)
        ]
        if finite:
            mean_metrics[str(metric)] = float(f"{np.mean(finite):.4}")

    # carry over only recognized run metadata (it describes the source run, not the
    # computed metrics, so it stays valid on recalculation)
    existing = model.metrics.get("diatomics", {})
    existing = existing if isinstance(existing, dict) else {}
    run_metadata = run_metadata or {}
    block: dict[str, DiatomicsYamlValue] = {}
    existing_pred_file = existing.get("pred_file")
    pred_file = existing_pred_file
    if pred_file_path is not None:
        pred_file = repo_relative_path(pred_file_path)
    if pred_file is not None:
        block["pred_file"] = pred_file

    pred_file_url = run_metadata.get("pred_file_url")
    if pred_file_url is None and pred_file == existing_pred_file:
        pred_file_url = existing.get("pred_file_url")
    if pred_file_url is not None:
        block["pred_file_url"] = pred_file_url

    for key in ("hardware", "run_time_sec", "excluded_formulas"):
        val = run_metadata.get(key)
        if val is None:
            val = existing.get(key)
        if val is not None:
            block[key] = val
    block |= mean_metrics

    # preserve_existing=False so a recompute fully replaces the block, dropping
    # deprecated metrics left in the YAML. This still runs when no finite metrics were
    # produced, preventing stale metric values from surviving.
    update_yaml_file(
        model.yaml_path, "metrics.diatomics", block, preserve_existing=False
    )
    print(f"Wrote {model.label} diatomic metrics to {model.yaml_path}")
    return block
