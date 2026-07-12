"""Metrics for evaluating ML thermal conductivity against DFT references.

Included metrics:
- SRD (Symmetric Relative Difference): Relative difference between ML and DFT
  conductivities
- SRE (Symmetric Relative Error): Absolute value of SRD
- SRME (Symmetric Relative Mean Error): Microscopic, mode-resolved error that avoids
  cancellation between over- and underpredicted mode contributions.

Adapted from https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/benchmark.py,
published in https://arxiv.org/abs/2408.00755. Ported in
https://github.com/janosh/matbench-discovery/pull/196 to parallelize over input
structures and scale thermal-conductivity metrics to larger test sets.
"""

import traceback
import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.enums import DataFiles, MbdKey, Model
from matbench_discovery.phonons import thermal_conductivity as ltc

KAPPA_ERROR_MAX = 2.0


class InvalidKappaReferenceError(ValueError):
    """Raised when reference conductivity data cannot be scored safely."""


def _finite_array(value: object, *, nonnegative: bool = False) -> np.ndarray | None:
    """Return finite numeric data, or None for malformed values."""
    try:
        values = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return None
    if values.size == 0 or np.any(~np.isfinite(values)):
        return None
    if nonnegative and np.any(values < 0):
        return None
    return values


def _single_kappa(value: object, *, reference: bool = False) -> float | None:
    """Return one non-negative conductivity, failing hard for bad references."""
    values = _finite_array(value, nonnegative=True)
    if values is None or values.size != 1:
        if reference:
            raise InvalidKappaReferenceError(
                f"Reference kappa must be one finite non-negative value, got {value!r}"
            )
        return None
    return float(values.ravel()[0])


def _symmetric_relative_difference(predicted: object, reference: object) -> float:
    """Return bounded SRD, assigning maximum error to invalid predictions."""
    true_value = _single_kappa(reference, reference=True)
    if true_value is None:  # pragma: no cover - reference=True raises instead
        raise ValueError("Invalid reference kappa")
    pred_value = _single_kappa(predicted)
    if pred_value is None:
        return -KAPPA_ERROR_MAX
    denominator = pred_value + true_value
    return (
        0.0
        if denominator == 0
        else KAPPA_ERROR_MAX * (pred_value - true_value) / denominator
    )


def calc_kappa_metrics_from_dfs(
    df_pred: pd.DataFrame,
    df_true: pd.DataFrame,
    *,
    ignore_imaginary_freqs: bool = False,
) -> pd.DataFrame:
    """Compute per-material thermal-conductivity metrics from two dataframes.

    Adds averaged conductivities and SRD/SRE/SRME metrics to raw ML predictions,
    handling array-valued columns such as tensors and mode-resolved properties.

    Args:
        df_pred: ML predictions with conductivity tensors, mode-resolved properties,
            and structural information.
        df_true: DFT references with the same structure as df_pred.
        ignore_imaginary_freqs: Score finite conductivity results even when their
            frequency diagnostic reports imaginary modes, matching protocols that
            explicitly opted into this historical benchmark behavior.

    Returns:
        df_pred with added benchmark columns:
        - SRD: Symmetric Relative Difference between ML and DFT conductivities
        - SRE: Absolute value of SRD
        - SRME: Mode-resolved error
        - DFT_kappa_tot_avg: Reference DFT conductivity values
    """
    df_pred = df_pred.copy()
    if ignore_imaginary_freqs:
        df_pred[Key.has_imag_ph_modes] = False

    df_pred[MbdKey.kappa_tot_avg] = df_pred[MbdKey.kappa_tot_rta].map(
        calculate_kappa_avg
    )

    true_kappa = df_true.loc[df_pred.index, MbdKey.kappa_tot_avg]
    df_pred[Key.srd] = [
        _symmetric_relative_difference(predicted, reference)
        for predicted, reference in zip(
            df_pred[MbdKey.kappa_tot_avg], true_kappa, strict=True
        )
    ]
    df_pred[Key.sre] = df_pred[Key.srd].abs()
    df_pred[Key.srme] = calc_kappa_srme_dataframes(df_pred, df_true)
    df_pred[MbdKey.true_kappa_tot_avg] = true_kappa

    return df_pred


def evaluate_kappa_predictions(
    df_predictions: pd.DataFrame,
    *,
    ignore_imaginary_freqs: bool = False,
) -> dict[str, float]:
    """Evaluate complete PhononDB predictions and return aggregate metrics."""
    from matbench_discovery.phonons import read_kappa_json

    df_reference = read_kappa_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path)
    if not df_predictions.index.is_unique:
        raise ValueError("Cannot evaluate kappa predictions: duplicate material IDs")
    if set(df_predictions.index) != set(df_reference.index):
        raise ValueError("Cannot evaluate kappa predictions: reference IDs differ")
    df_metrics = calc_kappa_metrics_from_dfs(
        df_predictions.loc[df_reference.index],
        df_reference,
        ignore_imaginary_freqs=ignore_imaginary_freqs,
    )
    metrics = {
        "srme": float(df_metrics[Key.srme].mean()),
        "sre": float(df_metrics[Key.sre].mean()),
    }
    if any(
        not np.isfinite(value) or not 0 <= value <= KAPPA_ERROR_MAX
        for value in metrics.values()
    ):
        raise ValueError(f"Invalid aggregate kappa metrics: {metrics}")
    return metrics


def weighted_quantiles(
    values: np.ndarray, weights: np.ndarray | None, quantile_levels: np.ndarray
) -> np.ndarray:
    """Weighted quantiles of a 1D sample via the empirical inverse CDF.

    Uses numpy's ``method="inverted_cdf"`` definition: Q(u) is the smallest x with
    F(x) >= u. Integer weights are exactly equivalent to repeated samples, so spectra
    on weighted irreducible q-meshes (DFT) and expanded full meshes (ML) compare
    identically.

    Args:
        values: 1D sample values; NaNs are dropped with their weights.
        weights: Non-negative 1D weights matching values, or None for uniform weights.
        quantile_levels: Levels in [0, 1] at which to evaluate the inverse CDF.

    Returns:
        Quantile values, same length as quantile_levels.

    Raises:
        ValueError: If values is empty (or all-NaN), lengths mismatch, weights are
            non-finite/negative/zero-sum, or quantile levels are outside [0, 1].
    """
    values = np.ravel(np.asarray(values, dtype=float))
    weights = (
        np.ones_like(values)
        if weights is None
        else np.ravel(np.asarray(weights, dtype=float))
    )
    quantile_levels = np.ravel(np.asarray(quantile_levels, dtype=float))
    if len(values) != len(weights):
        raise ValueError(f"{len(values)=} != {len(weights)=}")
    if np.any(~np.isfinite(quantile_levels)) or np.any(
        (quantile_levels < 0) | (quantile_levels > 1)
    ):
        raise ValueError("quantile_levels must be finite values in [0, 1]")
    finite_mask = np.isfinite(values)
    values, weights = values[finite_mask], weights[finite_mask]
    if len(values) == 0:
        raise ValueError("no finite values to compute quantiles from")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite non-negative values")
    weight_sum = weights.sum()
    if weight_sum <= 0:
        raise ValueError("weights must sum to a positive value")
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights) / weight_sum
    idx = np.searchsorted(cdf, quantile_levels, side="left")
    return values[np.minimum(idx, len(values) - 1)]


def mode_weights_for_freqs(
    ph_freqs: np.ndarray, q_weights: np.ndarray | None
) -> np.ndarray | None:
    """Per-frequency weights for a (n_qpoints, n_modes) phonon frequency array.

    q-point weights apply only when their length matches the q-point axis (irreducible
    mesh). Otherwise, a uniform full mesh is assumed and None is returned.
    """
    if q_weights is None:
        return None
    q_weights = np.ravel(np.asarray(q_weights, dtype=float))
    if q_weights.shape != (ph_freqs.shape[0],):
        return None
    return np.repeat(q_weights, ph_freqs.shape[1])


def calculate_kappa_avg(kappa: np.ndarray) -> np.ndarray | float:
    """Directionally average thermal conductivity from WTE-RTA outputs.

    Returns the trace average for tensors, i.e. the mean of diagonal components across
    the 3 spatial directions, as a scalar metric for comparing materials.

    Args:
        kappa: Thermal conductivity tensor of shape (..., 3, 3), a Voigt 6-vector
            [xx, yy, zz, yz, xz, xy] of shape (..., 6) as stored by phono3py RTA in
            calc_kappa.py, or pre-averaged directional conductivities of shape (..., 3).
            Earlier dimensions may encode temperatures or other parameters.

    Returns:
        Average conductivity: scalar for a single 3x3 tensor, array for multiple
        temperatures, or np.array([np.nan]) if calculation fails.
    """
    try:
        kappa_arr = np.asarray(kappa, dtype=float)
        if kappa_arr.ndim == 0:  # scalar NaN from a failed calculation
            raise ValueError(f"expected array-like kappa, got {kappa_arr=}")
        if kappa_arr.shape[-2:] == (3, 3):
            return np.trace(kappa_arr, axis1=-2, axis2=-1) / 3
        if kappa_arr.shape[-1] == 6:  # Voigt: diagonal components come first
            return kappa_arr[..., :3].mean(axis=-1)
        if kappa_arr.shape[-1:] == (3,):
            return kappa_arr.mean(axis=-1)
        raise ValueError(
            f"expected shape (..., 3, 3), (..., 6) or (..., 3), got {kappa_arr.shape}"
        )
    except (ValueError, TypeError):
        warnings.warn(
            f"Failed to calculate kappa_avg: {traceback.format_exc()}", stacklevel=2
        )
        return np.array([np.nan])


def calc_kappa_srme_dataframes(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> list[float]:
    """Calculate Symmetric Relative Mean Error (SRME) for each material.

    SRME measures total thermal-conductivity accuracy and individual phonon-mode
    contributions. Like SRD, it is symmetric between over- and underprediction, and it
    averages errors over all modes weighted by their contributions.

    Edge cases:
    - Returns 2.0 for materials with imaginary frequencies (unphysical predictions)
    - Returns 2.0 for materials where symmetry is broken during relaxation
    - Returns 2.0 for failed calculations or missing data

    Args:
        df_pred: ML predictions with mode-resolved properties.
        df_true: DFT references with mode-resolved properties.

    Returns:
        SRME values for each material, between 0 and 2:
        - 0 indicates perfect agreement in both total κ and mode-resolved properties
        - 2 indicates complete failure (imaginary frequencies, broken symmetry, etc.)
        - Intermediate values indicate partial agreement, lower being better
    """
    srme_list: list[float] = []
    for row_idx, row_pred in df_pred.iterrows():
        row_true = df_true.loc[row_idx]

        if row_pred.get(Key.has_imag_ph_modes) is True:
            srme_list.append(KAPPA_ERROR_MAX)
            continue
        if relaxed_space_group_number := row_pred.get(Key.final_spg_num):
            initial_space_group_number = row_pred.get(Key.init_spg_num)
            if (
                initial_space_group_number
                and relaxed_space_group_number != initial_space_group_number
            ) or (
                not initial_space_group_number
                and relaxed_space_group_number != row_true.get(Key.spg_num)
            ):
                srme_list.append(KAPPA_ERROR_MAX)
                continue
        _single_kappa(row_true[MbdKey.kappa_tot_avg], reference=True)
        try:
            result = np.ravel(calc_kappa_srme(row_pred, row_true))
        except InvalidKappaReferenceError:
            raise
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            result = np.array([KAPPA_ERROR_MAX])
        score = float(result[0]) if result.size == 1 else KAPPA_ERROR_MAX
        srme_list.append(
            score
            if np.isfinite(score) and 0 <= score <= KAPPA_ERROR_MAX
            else KAPPA_ERROR_MAX
        )

    return srme_list


def calc_kappa_srme(kappas_pred: pd.Series, kappas_true: pd.Series) -> np.ndarray:
    """Calculate Symmetric Relative Mean Error (SRME) for one material.

        SRME = 2 * (sum|κ_pred,i - κ_true,i| * w_i) / (κ_pred,tot + κ_true,tot)

    where:
    - κ_pred,i and κ_true,i are mode-resolved conductivities for mode i
    - w_i are the mode weights
    - κ_pred,tot and κ_true,tot are total conductivities

    Steps:
    1. Compute mode-resolved average conductivities if not precomputed
    2. Compute weighted mean absolute error over all modes
    3. Normalize by total-conductivity sum, making the metric symmetric and relative

    Args:
        kappas_pred: ML predictions including:
            - kappa_tot_avg: Average total conductivity
            - mode_kappa_tot: Mode-resolved conductivities
            - mode_weights: Mode weights for averaging
        kappas_true: DFT references with the same structure

    Returns:
        SRME values per temperature, each between 0 and 2:
        - 0 indicates perfect agreement in both total κ and mode-resolved properties
        - 2 indicates complete disagreement or invalid results
        Missing data or NaNs return np.array([2.0]).
    """
    true_totals = _finite_array(kappas_true[MbdKey.kappa_tot_avg], nonnegative=True)
    if true_totals is None:
        raise InvalidKappaReferenceError(
            "Reference kappa totals must be finite and non-negative"
        )
    true_totals = np.atleast_1d(true_totals)
    failure = np.full(true_totals.shape, KAPPA_ERROR_MAX)
    mode_kappa_tot_avgs: dict[str, np.ndarray] = {}
    # Try different data sources in order of preference for both pred and true data
    for label, kappas in {"true": kappas_true, "preds": kappas_pred}.items():
        keys = set(kappas.index)
        if MbdKey.mode_kappa_tot_avg in keys:
            mode_kappa = kappas[MbdKey.mode_kappa_tot_avg]
        elif MbdKey.mode_kappa_tot_rta in keys:
            mode_kappa = calculate_kappa_avg(kappas[MbdKey.mode_kappa_tot_rta])
        elif {MbdKey.kappa_p_rta, MbdKey.kappa_c, Key.heat_capacity} <= keys:
            try:
                mode_kappa = calculate_kappa_avg(
                    ltc.calc_mode_kappa_tot(
                        kappas[MbdKey.kappa_p_rta],
                        kappas[MbdKey.kappa_c],
                        kappas[Key.heat_capacity],
                    )
                )
            except (TypeError, ValueError) as exc:
                if label == "true":
                    raise InvalidKappaReferenceError(
                        "Invalid reference mode conductivity components"
                    ) from exc
                mode_kappa = np.nan
        else:
            if label == "preds":
                return failure
            raise InvalidKappaReferenceError(
                f"Neither mode_kappa_tot_avg, mode_kappa_tot nor individual kappa\n"
                f"components found in {label}, got\n{keys}"
            )
        mode_values = _finite_array(mode_kappa)
        if mode_values is None:
            if label == "true":
                raise InvalidKappaReferenceError(
                    "Reference mode conductivities must be finite"
                )
            return failure
        mode_kappa_tot_avgs[label] = mode_values

    predicted_totals = _finite_array(
        kappas_pred.get(MbdKey.kappa_tot_avg), nonnegative=True
    )
    mode_weights = _finite_array(kappas_pred.get(Key.mode_weights), nonnegative=True)
    if (
        predicted_totals is None
        or np.atleast_1d(predicted_totals).shape != true_totals.shape
        or mode_weights is None
        or mode_weights.sum() <= 0
    ):
        return failure
    predicted_totals = np.atleast_1d(predicted_totals)
    if mode_kappa_tot_avgs["preds"].shape != mode_kappa_tot_avgs["true"].shape:
        return failure

    # calculating microscopic error for all temperatures
    microscopic_error = np.atleast_1d(
        np.abs(mode_kappa_tot_avgs["preds"] - mode_kappa_tot_avgs["true"]).sum(
            axis=tuple(range(1, mode_kappa_tot_avgs["preds"].ndim))
        )
        / mode_weights.sum()
    )
    if microscopic_error.shape != true_totals.shape:
        return failure
    denominator = predicted_totals + true_totals
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = np.where(
            denominator == 0,
            np.where(microscopic_error == 0, 0, KAPPA_ERROR_MAX),
            KAPPA_ERROR_MAX * microscopic_error / denominator,
        )
    valid = np.isfinite(scores) & (scores >= 0) & (scores <= KAPPA_ERROR_MAX)
    return np.where(valid, scores, KAPPA_ERROR_MAX)


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, float],
    pred_file_path: str,
    *,
    run_metadata: Mapping[str, object] | None = None,
    force_file_path: str | None = None,
    run_info_path: str | None = None,
    replace_pred_file: bool = False,
) -> None:
    """Write kappa metrics to a model YAML file's phonons section.

    Args:
        model: Model to write metrics for.
        metrics: Kappa metrics for this model.
        pred_file_path: Path to prediction file.
        run_metadata: Optional complete-run timing, hardware, and memory provenance.
        force_file_path: Optional separate force-set artifact path.
        run_info_path: Optional small manifest/provenance sidecar path.
        replace_pred_file: Replace an existing local prediction path and clear its
            stale remote URL. The legacy default preserves established artifact paths.
    """
    from matbench_discovery import repo_relative_path
    from matbench_discovery.data import update_yaml_file

    pred_file_path = repo_relative_path(pred_file_path)

    def update_kappa_103(kappa_103: dict[str, Any]) -> dict[str, Any]:
        """Update metrics and provenance from the latest locked YAML section."""
        kappa_103.update(
            κ_SRME=round(metrics["srme"], 4),
            κ_SRE=round(metrics["sre"], 4),
        )
        if replace_pred_file:
            kappa_103 |= {"pred_file": pred_file_path, "pred_file_url": None}
        else:
            kappa_103.setdefault("pred_file", pred_file_path)
        for artifact_key, artifact_path in (
            ("force_file", force_file_path),
            ("run_info_file", run_info_path),
        ):
            url_key = f"{artifact_key}_url"
            if artifact_path is not None:
                relative_path = repo_relative_path(artifact_path)
                artifact_changed = kappa_103.get(artifact_key) != relative_path
                kappa_103[artifact_key] = relative_path
                if replace_pred_file or artifact_changed:
                    kappa_103[url_key] = None
            elif replace_pred_file:
                kappa_103.pop(artifact_key, None)
                kappa_103.pop(url_key, None)
        for key in ("hardware", "run_time_sec", "max_rss_gb", "max_gpu_mem_gb"):
            if run_metadata is not None and key in run_metadata:
                kappa_103[key] = run_metadata[key]
            elif replace_pred_file:
                kappa_103.pop(key, None)
        return kappa_103

    update_yaml_file(
        model.yaml_path,
        "metrics.phonons.kappa_103",
        update_kappa_103,
        preserve_existing=False,
    )
