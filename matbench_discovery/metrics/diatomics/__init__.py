"""Metrics for analyzing diatomic potential energy curves.

Big props to Tamas Stenczel and Yuan Chiang who spearheaded this type of PES
smoothness analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP
paper https://arxiv.org/abs/2401.00096 (see fig. 56) and MLIP Arena
https://huggingface.co/spaces/atomind/mlip-arena, respectively.
"""

from dataclasses import dataclass
from typing import Any, Self

import numpy as np

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.diatomics import energy, force  # noqa: F401
from matbench_discovery.metrics.diatomics.energy import (
    calc_curve_diff_auc,
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae,
    calc_second_deriv_smoothness,
    calc_tortuosity,
)
from matbench_discovery.metrics.diatomics.force import (
    calc_conservation_deviation,
    calc_force_flips,
    calc_force_jump,
    calc_force_mae,
    calc_force_total_variation,
)


@dataclass
class DiatomicCurve:
    """Energies and forces for a single diatomic molecule at multiple distances."""

    distances: np.ndarray  # shape (n_distances,)
    energies: np.ndarray  # shape (n_distances,)
    forces: np.ndarray  # shape (n_distances, n_atoms, 3)

    def __post_init__(self) -> None:
        """Convert inputs to numpy arrays."""
        self.energies = np.asarray(self.energies)
        self.forces = np.asarray(self.forces)
        self.distances = np.asarray(self.distances)


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
    hetero_nuclear: dict[str, DiatomicCurve] | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create DiatomicCurves from a dictionary loaded from JSON."""
        dists = data["distances"] = np.asarray(data["distances"])
        for key in {"homo-nuclear", "hetero-nuclear"} & set(data):
            data[key.replace("-", "_")] = {
                k: DiatomicCurve(**v, distances=dists) for k, v in data.pop(key).items()
            }
        return cls(**data)


def calc_diatomic_curve_metrics(
    ref_curves: DiatomicCurves | None,
    pred_curves: DiatomicCurves,
    metrics: dict[str, dict[str, Any]] | None = None,
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

    Returns:
        dict[str, dict[str, float]]: Map of element symbols to metric dicts with keys
            being the metric names and values being the metric values.
    """
    if unknown_metrics := set(metrics or {}) - set(MbdKey):
        raise ValueError(f"{unknown_metrics=}. Valid metrics=")

    results: dict[str, dict[str, float]] = {}
    metrics = (metrics or {}).copy()

    # Initialize empty kwargs for each metric if not provided
    for key in MbdKey:
        metrics.setdefault(key, {})

    for elem_symbol, pred_data in pred_curves.homo_nuclear.items():
        elem_metrics: dict[str, float] = {}
        distances = pred_curves.distances

        # Skip reference-requiring metrics if no reference curves provided
        if ref_curves and (ref_data := ref_curves.homo_nuclear.get(elem_symbol)):
            if not np.array_equal(distances, ref_curves.distances):
                raise ValueError(
                    "Reference and predicted distances must be the same. If goal is "
                    "to interpolate predicted curves to reference distances, do so "
                    "before passing to calc_diatomic_curve_metrics."
                )

            # Energy metrics that need both curves
            if MbdKey.norm_auc in metrics:
                elem_metrics[MbdKey.norm_auc] = calc_curve_diff_auc(
                    distances,
                    ref_data.energies,
                    distances,
                    pred_data.energies,
                    **metrics[MbdKey.norm_auc],
                )

            if MbdKey.energy_mae in metrics:
                elem_metrics[MbdKey.energy_mae] = calc_energy_mae(
                    distances,
                    ref_data.energies,
                    distances,
                    pred_data.energies,
                    **metrics[MbdKey.energy_mae],
                )

            if MbdKey.force_mae in metrics:
                elem_metrics[MbdKey.force_mae] = calc_force_mae(
                    distances,
                    ref_data.forces,
                    distances,
                    pred_data.forces,
                    **metrics[MbdKey.force_mae],
                )

        # Energy metrics that need only predicted curve
        if MbdKey.smoothness in metrics:
            elem_metrics[MbdKey.smoothness] = calc_second_deriv_smoothness(
                distances, pred_data.energies, **metrics[MbdKey.smoothness]
            )

        if MbdKey.tortuosity in metrics:
            elem_metrics[MbdKey.tortuosity] = calc_tortuosity(
                distances, pred_data.energies, **metrics[MbdKey.tortuosity]
            )

        if MbdKey.energy_diff_flips in metrics:
            elem_metrics[MbdKey.energy_diff_flips] = calc_energy_diff_flips(
                distances, pred_data.energies, **metrics[MbdKey.energy_diff_flips]
            )

        if MbdKey.energy_grad_norm_max in metrics:
            elem_metrics[MbdKey.energy_grad_norm_max] = calc_energy_grad_norm_max(
                distances, pred_data.energies, **metrics[MbdKey.energy_grad_norm_max]
            )

        if MbdKey.energy_jump in metrics:
            elem_metrics[MbdKey.energy_jump] = calc_energy_jump(
                distances, pred_data.energies, **metrics[MbdKey.energy_jump]
            )

        # Force metrics that need only predicted curve
        if MbdKey.conservation in metrics:
            elem_metrics[MbdKey.conservation] = calc_conservation_deviation(
                distances,
                pred_data.energies,
                pred_data.forces,
                **metrics[MbdKey.conservation],
            )

        if MbdKey.force_flips in metrics:
            elem_metrics[MbdKey.force_flips] = calc_force_flips(
                distances, pred_data.forces, **metrics[MbdKey.force_flips]
            )

        if MbdKey.force_total_variation in metrics:
            elem_metrics[MbdKey.force_total_variation] = calc_force_total_variation(
                distances, pred_data.forces, **metrics[MbdKey.force_total_variation]
            )

        if MbdKey.force_jump in metrics:
            elem_metrics[MbdKey.force_jump] = calc_force_jump(
                distances, pred_data.forces, **metrics[MbdKey.force_jump]
            )

        results[elem_symbol] = elem_metrics

    return results
