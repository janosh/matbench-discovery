"""Metrics for analyzing diatomic potential energy curves.

Big props to Tamas Stenczel and Yuan Chiang who spearheaded this type of PES
smoothness analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP
paper https://arxiv.org/abs/2401.00096 (see fig. 56) and MLIP Arena
https://huggingface.co/spaces/atomind/mlip-arena, respectively.
"""

import inspect
from collections.abc import Callable, Mapping, Sequence
from typing import Any

# ruff: noqa: F401 (don't flag convenience imports above)
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.diatomics import energy, force
from matbench_discovery.metrics.diatomics.energy import (
    calc_conservation_deviation,
    calc_curve_diff_auc,
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae_vs_ref,
    calc_second_deriv_smoothness,
    calc_tortuosity,
)
from matbench_discovery.metrics.diatomics.force import (
    calc_force_flips,
    calc_force_jump,
    calc_force_mae_vs_ref,
    calc_force_total_variation,
)

# Type alias for a curve represented as a tuple of x and y values
DiatomicCurve = tuple[Sequence[float], Sequence[float]]
# Type alias for a dictionary mapping element symbols to curves
DiatomicCurves = Mapping[str, DiatomicCurve]


def calc_diatomic_curve_metrics(
    ref_curves: DiatomicCurves,
    pred_curves: DiatomicCurves,
    pred_force_curves: DiatomicCurves | None = None,
    metrics: dict[str, dict[str, Any]] | None = None,
) -> dict[str, dict[str, float]]:
    """Calculate diatomic curve metrics comparing predicted curves to reference curves.

    Args:
        ref_curves (DiatomicCurves): Reference energy curves for each element.
        pred_curves (DiatomicCurves): Predicted energy curves for each element.
        pred_force_curves (DiatomicCurves | None): Predicted force curves for each
            element. Required for force-based metrics.
        metrics (dict[str, dict[str, Any]] | None): Dictionary mapping metric names to
            dictionaries of keyword arguments for each metric function. If None, uses
            all metrics with default parameters. To use a subset of metrics, provide
            a dictionary with those metric names as keys and their keyword arguments
            as values. Empty dictionaries will use default parameters.

    Returns:
        dict[str, dict[str, float]]: Dictionary mapping element symbols to metrics dict
            with keys being the metric names and values being the metric values.
    """
    results: dict[str, dict[str, float]] = {}

    # Map metric keys to their functions
    metric_functions: dict[str, Callable[..., float]] = {
        # Energy metrics that need both curves
        MbdKey.norm_auc: energy.calc_curve_diff_auc,
        MbdKey.energy_mae_vs_ref: energy.calc_energy_mae_vs_ref,
        # Energy metrics that need only predicted curve
        MbdKey.smoothness: energy.calc_second_deriv_smoothness,
        MbdKey.tortuosity: energy.calc_tortuosity,
        MbdKey.conservation: energy.calc_conservation_deviation,
        MbdKey.energy_diff_flips: energy.calc_energy_diff_flips,
        MbdKey.energy_grad_norm_max: energy.calc_energy_grad_norm_max,
        MbdKey.energy_jump: energy.calc_energy_jump,
        # Force metrics that need both curves
        MbdKey.force_mae_vs_ref: force.calc_force_mae_vs_ref,
        # Force metrics that need only predicted curve
        MbdKey.force_flips: force.calc_force_flips,
        MbdKey.force_total_variation: force.calc_force_total_variation,
        MbdKey.force_jump: force.calc_force_jump,
    }

    # If no metrics specified, use all metrics with default parameters
    metrics = (metrics or {}).copy()

    if unknown_metrics := set(metrics) - set(metric_functions):
        raise ValueError(
            f"{unknown_metrics=}. Valid metrics={', '.join(metric_functions)}"
        )

    for key in metric_functions:
        metrics.setdefault(key, {})

    # Remove force-based metrics if no force curves provided
    if pred_force_curves is None:
        metrics = {
            name: kwargs
            for name, kwargs in metrics.items()
            if not name.startswith("force_")
        }

    for elem_symbol, ref_curve in ref_curves.items():
        if elem_symbol not in pred_curves:
            continue

        pred_curve = pred_curves[elem_symbol]
        elem_metrics: dict[str, float] = {}

        for name, func_kwargs in metrics.items():
            metric_func = metric_functions[name]
            param_set = set(inspect.signature(metric_func).parameters)
            needs_ref_curve = any("_ref" in param for param in param_set)
            is_force_metric = name.startswith("force_")

            # Handle force metrics
            if is_force_metric:
                if pred_force_curves is None or elem_symbol not in pred_force_curves:
                    continue
                curve_to_use = pred_force_curves[elem_symbol]
            else:  # energy metrics
                curve_to_use = pred_curve

            # Call metric function with appropriate arguments
            if needs_ref_curve:
                elem_metrics[name] = metric_func(
                    *ref_curve, *curve_to_use, **func_kwargs
                )
            else:
                elem_metrics[name] = metric_func(*curve_to_use, **func_kwargs)

        results[elem_symbol] = elem_metrics

    return results
