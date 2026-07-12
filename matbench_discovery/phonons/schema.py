"""Canonicalize thermal-conductivity records from current and legacy runners."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey

if TYPE_CHECKING:
    from collections.abc import Mapping

KAPPA_TENSOR_KEYS = (
    str(MbdKey.kappa_tot_rta),
    str(MbdKey.kappa_p_rta),
    str(MbdKey.kappa_c),
    str(MbdKey.mode_kappa_tot_rta),
)
VOIGT_MATRIX_INDICES = np.array(((0, 5, 4), (5, 1, 3), (4, 3, 2)))

KAPPA_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    str(Key.mat_id): ("mp_id",),
    str(Key.init_spg_num): (
        str(Key.spg_num),
        "initial_space_group_number",
        "initial_spg_num",
        "init_space_group_number",
    ),
    str(Key.final_spg_num): (
        "relaxed_space_group_number",
        "final_space_group_number",
        "relaxed_spg_num",
    ),
    str(Key.has_imag_ph_modes): (
        "imaginary_freqs",
        "has_imaginary_freqs",
        "has_imaginary_modes",
    ),
    str(Key.ph_freqs): ("frequencies", "phonon_frequencies"),
}


def _is_missing_scalar(value: object) -> bool:
    """Return whether a scalar is a DataFrame placeholder rather than alias data."""
    missing = pd.isna(value)
    return isinstance(missing, bool | np.bool_) and bool(missing)


def _has_errors(value: object) -> bool:
    """Treat missing error cells as empty and populated values as failures."""
    return not _is_missing_scalar(value) and bool(value)


def voigt_6_to_full_3x3(tensor: object) -> object:
    """Expand a trailing Voigt-6 axis to a symmetric 3x3 tensor.

    Arrays whose trailing dimension is not six are returned unchanged.
    """
    try:
        tensor_array = np.asarray(tensor, dtype=float)
    except (TypeError, ValueError):
        return tensor
    if tensor_array.ndim == 0 or tensor_array.shape[-1] != 6:
        return tensor
    return tensor_array[..., VOIGT_MATRIX_INDICES]


def _values_equal(first: object, second: object) -> bool:
    """Compare scalar or array alias values, treating paired NaNs as equal."""
    try:
        first_array, second_array = np.asarray(first), np.asarray(second)
        if first_array.shape != second_array.shape:
            return False
        try:
            return bool(np.array_equal(first_array, second_array, equal_nan=True))
        except TypeError:
            return bool(np.array_equal(first_array, second_array))
    except (TypeError, ValueError):
        return False


def normalize_kappa_result(result: Mapping[str, Any]) -> dict[str, Any]:
    """Return one result row using canonical IDs, symmetry fields, and tensors."""
    normalized = {str(key): value for key, value in result.items()}
    for canonical, aliases in KAPPA_COLUMN_ALIASES.items():
        present_names = [
            name
            for name in (canonical, *aliases)
            if name in normalized and not _is_missing_scalar(normalized[name])
        ]
        if present_names:
            value = normalized[present_names[0]]
            if any(
                not _values_equal(value, normalized[name]) for name in present_names[1:]
            ):
                raise ValueError(f"Conflicting aliases for kappa column {canonical!r}")
            normalized[canonical] = value
        for alias in aliases:
            normalized.pop(alias, None)

    for tensor_key in KAPPA_TENSOR_KEYS:
        if tensor_key in normalized:
            normalized[tensor_key] = voigt_6_to_full_3x3(normalized[tensor_key])

    initial_spg = normalized.get(str(Key.init_spg_num))
    final_spg = normalized.get(str(Key.final_spg_num))
    if "broken_symmetry" not in normalized:
        normalized["broken_symmetry"] = bool(
            not pd.isna(initial_spg)
            and not pd.isna(final_spg)
            and initial_spg != final_spg
        )
    normalized.setdefault("errors", [])
    normalized.setdefault("error_traceback", [])
    if _has_errors(normalized["errors"]) or normalized.get("conductivity_skipped"):
        normalized[str(MbdKey.kappa_tot_rta)] = np.nan
        normalized[str(MbdKey.mode_kappa_tot_rta)] = np.nan
        normalized[str(Key.mode_weights)] = np.nan
    return normalized


def normalize_kappa_dataframe(df_kappa: pd.DataFrame) -> pd.DataFrame:
    """Normalize legacy kappa dataframe columns without changing row order."""
    normalized = df_kappa.copy()
    if str(Key.mat_id) not in normalized and normalized.index.name in {
        str(Key.mat_id),
        *KAPPA_COLUMN_ALIASES[str(Key.mat_id)],
    }:
        normalized = normalized.reset_index()
    rows = [normalize_kappa_result(row) for row in normalized.to_dict(orient="records")]
    return pd.DataFrame(rows, index=normalized.index)
