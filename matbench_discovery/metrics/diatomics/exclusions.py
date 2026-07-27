"""Model-specific exclusions for diatomic metric aggregation."""

from matbench_discovery.enums import Model
from matbench_discovery.metrics.diatomics import NON_MP_ELEMENTS


def is_non_mp_formula(formula: str) -> bool:
    """Whether a diatomic formula involves an element outside the MP element set."""
    return any(element in NON_MP_ELEMENTS for element in formula.split("-"))


def get_excluded_formula_reasons(
    model_key: str, invalid_formulas: tuple[str, ...] | list[str] = ()
) -> dict[str, str]:
    """Combine YAML-curated and run-discovered exclusions for one model.

    Curated reasons take precedence. Formulas containing elements outside the
    Materials Project element set are omitted because metrics skip them globally.
    """
    try:
        diatomics_metrics = Model.from_ref(model_key).metrics.get("diatomics") or {}
    except ValueError:  # debug models like emt have no Model enum entry
        diatomics_metrics = {}
    curated_reasons = diatomics_metrics.get("excluded_formula_reasons", {})
    reasons = dict.fromkeys(invalid_formulas, "invalid or unsupported curve")
    reasons |= curated_reasons
    return {
        formula: reasons[formula]
        for formula in sorted(reasons)
        if not is_non_mp_formula(formula)
    }


def drop_metric_exclusions(
    model_key: str, metrics: dict[str, dict[str, float]]
) -> dict[str, dict[str, float]]:
    """Remove model-specific pathological curves before metric aggregation."""
    excluded = set(get_excluded_formula_reasons(model_key))
    return {
        key: value
        for key, value in metrics.items()
        if key not in excluded and f"{key}-{key}" not in excluded
    }
