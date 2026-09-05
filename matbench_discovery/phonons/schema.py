"""Validate canonical thermal-conductivity records and normalize internal tensors."""

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

# Reject retired spellings instead of silently ignoring data in legacy fields.
NON_CANONICAL_FIELDS = frozenset(
    {
        "mp_id",
        "spg_num",
        "initial_space_group_number",
        "initial_spg_num",
        "init_space_group_number",
        "relaxed_space_group_number",
        "final_space_group_number",
        "relaxed_spg_num",
        "imaginary_freqs",
        "has_imaginary_freqs",
        "has_imaginary_modes",
        "frequencies",
        "phonon_frequencies",
    }
)


def _is_missing_scalar(value: object) -> bool:
    """Return whether a scalar is a missing DataFrame value."""
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
    except TypeError, ValueError:
        return tensor
    if tensor_array.ndim == 0 or tensor_array.shape[-1] != 6:
        return tensor
    return tensor_array[..., VOIGT_MATRIX_INDICES]


def normalize_kappa_result(result: Mapping[str, Any]) -> dict[str, Any]:
    """Return one result row using canonical IDs, symmetry fields, and tensors."""
    normalized = {str(key): value for key, value in result.items()}
    if legacy := normalized.keys() & NON_CANONICAL_FIELDS:
        raise ValueError(
            f"Non-canonical kappa fields {sorted(legacy)}; migrate the artifact"
        )
    material_id = normalized.get(str(Key.mat_id))
    if not isinstance(material_id, str) or not material_id.strip():
        raise ValueError(
            f"Kappa record requires a non-empty material_id, got {material_id!r}"
        )

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
    """Normalize canonical records without changing dataframe row order."""
    rows = [normalize_kappa_result(row) for row in df_kappa.to_dict(orient="records")]
    return pd.DataFrame(rows, index=df_kappa.index)
