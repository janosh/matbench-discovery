"""Metrics for evaluating ML thermal conductivities against DFT reference values.

The metrics include:
- SRD (Symmetric Relative Difference): Measures the relative difference between ML and
    DFT predictions
- SRE (Symmetric Relative Error): Absolute value of SRD
- SRME (Symmetric Relative Mean Error): A microscopic (i.e. mode-resolved)
  metric that is particularly useful since it's not subject to error cancellation.
  Overpredictions of kappa-contributions from one mode will not cancel against
  underpredictions from another mode.

Code in this module is adapted from https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/benchmark.py.
All credit to Balázs Póta, Paramvir Ahlawat, Gábor Csányi, Michele Simoncelli. See
https://arxiv.org/abs/2408.00755 for details.
It was ported to this repo in https://github.com/janosh/matbench-discovery/pull/196 to
implement parallelization across input structures which allows scaling thermal
conductivity metric to larger test sets.
"""

import traceback
import warnings

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.phonons import thermal_conductivity as ltc


def calc_kappa_metrics_from_dfs(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> pd.DataFrame:
    """Compute per-material thermal conductivity predictions metrics from 2 dataframes.

    This function takes the raw ML predictions and DFT reference results and computes
    various benchmark metrics. It handles array-type columns (like stress tensors and
    mode-resolved properties), calculates averaged quantities, and computes the
    SRD, SRE, and SRME metrics.

    Args:
        df_pred: DataFrame containing ML model predictions with columns for
            thermal conductivity tensors, mode-resolved properties, and other
            structural information.
        df_true: DataFrame containing DFT reference calculations with the
            same structure as df_pred.

    Returns:
        pd.DataFrame: df_pred with additional columns for benchmark metrics:
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


def calculate_kappa_avg(kappa: np.ndarray) -> np.ndarray:
    """Calculate directionally averaged trace of the conductivity tensor obtained from
    the Wigner transport equation (WTE) solution in the relaxation-time approximation.

    Takes a thermal conductivity tensor and returns its trace (average of diagonal
    components). This represents the average thermal conductivity in the 3 spatial
    directions, which is a useful scalar metric for comparing materials.

    Args:
        kappa: Thermal conductivity tensor, typically of shape (..., 3, 3) where
            the last two dimensions represent the 3x3 conductivity tensor.
            Earlier dimensions may include temperatures or other parameters.

    Returns:
        np.ndarray: Average conductivity value(s). Returns np.nan if the input contains
        any NaN values or if the calculation fails. For multiple temperatures,
        returns an array of averages.
    """
    if np.any(pd.isna(kappa)):
        return np.array([np.nan])
    try:
        return np.asarray(kappa)[..., :3].mean(axis=-1)
    except Exception:
        warnings.warn(
            f"Failed to calculate kappa_avg: {traceback.format_exc()}", stacklevel=2
        )
        return np.array([np.nan])


def calc_kappa_srme_dataframes(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> list[float]:
    """Calculate the Symmetric Relative Mean Error (SRME) for each material.

    SRME is a comprehensive metric that evaluates both the overall accuracy of thermal
    conductivity predictions and the accuracy of individual phonon mode contributions.
    It is symmetric (like SRD) to treat over- and under-predictions equally, and
    accounts for the mean error across all phonon modes weighted by their contributions.

    The function handles various edge cases:
    - Returns 2.0 for materials with imaginary frequencies (unphysical predictions)
    - Returns 2.0 for materials where symmetry is broken during relaxation
    - Returns 2.0 for failed calculations or missing data

    Args:
        df_pred (pd.DataFrame): ML predictions including mode-resolved properties
        df_true (pd.DataFrame): DFT reference data including mode-resolved properties

    Returns:
        list[float]: SRME values for each material. Values are between 0 and 2, where:
        - 0 indicates perfect agreement in both total κ and mode-resolved properties
        - 2 indicates complete failure (imaginary frequencies, broken symmetry, etc.)
        - Values in between indicate partial agreement, with lower being better
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
            if initial_space_group_number := row_pred.get(Key.init_spg_num):
                if relaxed_space_group_number != initial_space_group_number:
                    srme_list.append(2)
                    continue
            elif relaxed_space_group_number != row_true.get(Key.spg_num):
                srme_list.append(2)
                continue
        result = calc_kappa_srme(row_pred, row_true)
        srme_list.append(result[0])  # append the first temperature's SRME

    return srme_list


def calc_kappa_srme(kappas_pred: pd.Series, kappas_true: pd.Series) -> np.ndarray:
    """Calculate the Symmetric Relative Mean Error (SRME) for a single material.

        SRME = 2 * (sum|κ_pred,i - κ_true,i| * w_i) / (κ_pred,tot + κ_true,tot)

    where:
    - κ_pred,i and κ_true,i are mode-resolved conductivities for mode i
    - w_i are the mode weights
    - κ_pred,tot and κ_true,tot are total conductivities

    The calculation involves:
    1. Computing mode-resolved average conductivities if not pre-computed
    2. Calculating the weighted mean absolute error across all modes
    3. Normalizing by the sum of total conductivities to make it symmetric and relative

    Args:
        kappas_pred: Series containing ML predictions including:
            - kappa_tot_avg: Average total conductivity
            - mode_kappa_tot: Mode-resolved conductivities
            - mode_weights: Mode weights for averaging
        kappas_true: Series containing DFT reference data with same structure

    Returns:
        np.ndarray: SRME values per temperature, each between 0 and 2, where:
        - 0 indicates perfect agreement in both total κ and mode-resolved properties
        - 2 indicates complete disagreement or invalid results
        On error conditions (missing data, NaN values), returns np.array([2.0]).
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
            kappas = kappas[MbdKey.mode_kappa_tot_avg]
        elif MbdKey.mode_kappa_tot_rta in kappas:
            kappas = calculate_kappa_avg(kappas[MbdKey.mode_kappa_tot_rta])
        elif {MbdKey.kappa_p_rta, MbdKey.kappa_c, Key.heat_capacity} <= keys:
            kappas = calculate_kappa_avg(
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
        mode_kappa_tot_avgs[label] = np.asarray(kappas)

    # calculating microscopic error for all temperatures
    microscopic_error = (
        np.abs(mode_kappa_tot_avgs["preds"] - mode_kappa_tot_avgs["true"]).sum(
            axis=tuple(range(1, np.asarray(mode_kappa_tot_avgs["preds"]).ndim))
        )
        / np.asarray(kappas_pred[Key.mode_weights]).sum()
    )

    denominator = kappas_pred[MbdKey.kappa_tot_avg] + kappas_true[MbdKey.kappa_tot_avg]
    return 2 * microscopic_error / denominator
