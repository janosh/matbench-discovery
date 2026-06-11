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

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.phonons import thermal_conductivity as ltc


def calc_kappa_metrics_from_dfs(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> pd.DataFrame:
    """Compute per-material thermal-conductivity metrics from two dataframes.

    Adds averaged conductivities and SRD/SRE/SRME metrics to raw ML predictions,
    handling array-valued columns such as tensors and mode-resolved properties.

    Args:
        df_pred: ML predictions with conductivity tensors, mode-resolved properties,
            and structural information.
        df_true: DFT references with the same structure as df_pred.

    Returns:
        df_pred with added benchmark columns:
        - SRD: Symmetric Relative Difference between ML and DFT conductivities
        - SRE: Absolute value of SRD
        - SRME: Mode-resolved error
        - DFT_kappa_tot_avg: Reference DFT conductivity values
    """
    # Remove precomputed columns
    df_pred[MbdKey.kappa_tot_avg] = df_pred[MbdKey.kappa_tot_rta].map(
        calculate_kappa_avg
    )

    df_pred[Key.srd] = (
        2
        * (df_pred[MbdKey.kappa_tot_avg] - df_true[MbdKey.kappa_tot_avg])
        / (df_pred[MbdKey.kappa_tot_avg] + df_true[MbdKey.kappa_tot_avg])
    )

    # turn temperature list to the first temperature (300K) TODO: allow multiple
    # temperatures to be tested
    df_pred[Key.srd] = df_pred[Key.srd].map(
        lambda x: x if isinstance(x, float) else x[0]
    )

    # We substitute NaN values with 0 predicted conductivity, yielding -2 for SRD
    df_pred[Key.srd] = df_pred[Key.srd].fillna(-2)

    df_pred[Key.sre] = df_pred[Key.srd].abs()

    df_pred[Key.srme] = calc_kappa_srme_dataframes(df_pred, df_true)

    df_pred[MbdKey.true_kappa_tot_avg] = df_true[MbdKey.kappa_tot_avg]

    return df_pred


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
    for idx, row_pred in df_pred.iterrows():
        row_true = df_true.loc[idx]

        # NOTE code below just until before return used to be wrapped in try/except in
        # which case SRME=2 was set for the failing material
        if row_pred.get(Key.has_imag_ph_modes) is True:
            srme_list.append(2)
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
                srme_list.append(2)
                continue
        result = calc_kappa_srme(row_pred, row_true)
        srme_list.append(float(np.ravel(result)[0]))

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
    if np.any(np.isnan(kappas_true[MbdKey.kappa_tot_avg])):
        raise ValueError("found NaNs in kappa_tot_avg reference values")
    if (  # return highest possible SRME=2 if any of these conditions are met:
        # only have NaN averaged kappa preds
        np.all(np.isnan(kappas_pred[MbdKey.kappa_tot_avg]))
        # some mode-resolved kappa preds are NaN
        or np.any(np.isnan(kappas_pred[MbdKey.kappa_tot_rta]))
        # some mode weights are NaN
        or np.any(np.isnan(kappas_pred[Key.mode_weights]))
    ):
        return np.array([2.0])

    mode_kappa_tot_avgs = {}  # store results for pred and true
    # Try different data sources in order of preference for both pred and true data
    for label, kappas in {"preds": kappas_pred, "true": kappas_true}.items():
        keys = set(kappas.keys())
        if MbdKey.mode_kappa_tot_avg in kappas:
            mode_kappa = kappas[MbdKey.mode_kappa_tot_avg]
        elif MbdKey.mode_kappa_tot_rta in kappas:
            mode_kappa = calculate_kappa_avg(kappas[MbdKey.mode_kappa_tot_rta])
        elif {MbdKey.kappa_p_rta, MbdKey.kappa_c, Key.heat_capacity} <= keys:
            mode_kappa = calculate_kappa_avg(
                ltc.calc_mode_kappa_tot(
                    kappas[MbdKey.kappa_p_rta],
                    kappas[MbdKey.kappa_c],
                    kappas[Key.heat_capacity],
                )
            )
        else:
            raise ValueError(
                f"Neither mode_kappa_tot_avg, mode_kappa_tot nor individual kappa\n"
                f"components found in {label}, got\n{keys}"
            )
        mode_kappa_tot_avgs[label] = np.asarray(mode_kappa)

    # calculating microscopic error for all temperatures
    microscopic_error = (
        np.abs(mode_kappa_tot_avgs["preds"] - mode_kappa_tot_avgs["true"]).sum(
            axis=tuple(range(1, np.asarray(mode_kappa_tot_avgs["preds"]).ndim))
        )
        / np.asarray(kappas_pred[Key.mode_weights]).sum()
    )

    denominator = kappas_pred[MbdKey.kappa_tot_avg] + kappas_true[MbdKey.kappa_tot_avg]
    return 2 * microscopic_error / denominator


def write_metrics_to_yaml(
    model: Model, metrics: dict[str, float], pred_file_path: str
) -> None:
    """Write kappa metrics to a model YAML file's phonons section.

    Args:
        model: Model to write metrics for.
        metrics: Kappa metrics for this model.
        pred_file_path: Path to prediction file.
    """
    from matbench_discovery import ROOT
    from matbench_discovery import data as mbd_data

    # Convert absolute path to repo-relative path
    pred_file_path = pred_file_path.removeprefix(f"{ROOT}/")

    with open(model.yaml_path, encoding="utf-8") as file:
        data = mbd_data.round_trip_yaml.load(file)

    # Ensure nested structure exists and update non-destructively
    kappa_103 = (
        data.setdefault("metrics", {})
        .setdefault("phonons", {})
        .setdefault("kappa_103", {})
    )
    kappa_103.update(
        κ_SRME=float(round(metrics["srme"], 4)),
        κ_SRE=float(round(metrics["sre"], 4)),
    )
    kappa_103.setdefault("pred_file", pred_file_path)

    with open(model.yaml_path, mode="w", encoding="utf-8") as file:
        mbd_data.round_trip_yaml.dump(data, file)
