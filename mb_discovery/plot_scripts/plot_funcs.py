from __future__ import annotations

from typing import Any, Literal, Sequence, get_args

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.stats
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"

StabilityCriterion = Literal["energy", "energy+std", "energy-std"]

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=200)
plt.rc("font", size=16)


def hist_classified_stable_as_func_of_hull_dist(
    df: pd.DataFrame,
    target_col: str,
    pred_cols: Sequence[str],
    e_above_hull_col: str,
    ax: plt.Axes = None,
    energy_type: Literal["true", "pred"] = "true",
    stability_crit: StabilityCriterion = "energy",
    show_mae: bool = False,
    stability_threshold: float = 0,  # set stability threshold as distance to convex
    # hull in eV / atom, usually 0 or 0.1 eV
    x_lim: tuple[float, float] = (-0.4, 0.4),
) -> plt.Axes:
    """
    Histogram of the energy difference (either according to DFT ground truth [default]
    or model predicted energy) to the convex hull for materials in the WBM data set. The
    histogram is broken down into true positives, false negatives, false positives, and
    true negatives based on whether the model predicts candidates to be below the known
    convex hull. Ideally, in discovery setting a model should exhibit high recall, i.e.
    the majority of materials below the convex hull being correctly identified by the
    model.

    See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.

    NOTE this figure plots hist bars separately which causes aliasing in pdf
    to resolve this take into Inkscape and merge regions by color
    """
    if ax is None:
        ax = plt.gca()

    error = df[pred_cols].mean(axis=1) - df[target_col]
    e_above_hull_vals = df[e_above_hull_col]
    residuals = error + e_above_hull_vals

    if stability_crit not in get_args(StabilityCriterion):
        raise ValueError(
            f"Invalid {stability_crit=} must be one of {get_args(StabilityCriterion)}"
        )
    if stability_crit == "energy":
        test = residuals
    elif "std" in stability_crit:
        # TODO column names to compute standard deviation from are currently hardcoded
        # needs to be updated when adding non-aviary models with uncertainty estimation
        var_aleatoric = (df.filter(like="_ale_") ** 2).mean(axis=1)
        var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
        std_total = (var_epistemic + var_aleatoric) ** 0.5

        if stability_crit == "energy+std":
            test += std_total
        elif stability_crit == "energy-std":
            test -= std_total

    # --- histogram by DFT-computed distance to convex hull
    if energy_type == "true":
        actual_pos = e_above_hull_vals <= stability_threshold
        actual_neg = e_above_hull_vals > stability_threshold
        model_pos = test <= stability_threshold
        model_neg = test > stability_threshold

        n_true_pos = len(e_above_hull_vals[actual_pos & model_pos])
        n_false_neg = len(e_above_hull_vals[actual_pos & model_neg])

        n_total_pos = n_true_pos + n_false_neg
        null = n_total_pos / len(e_above_hull_vals)

        true_pos = e_above_hull_vals[actual_pos & model_pos]
        false_neg = e_above_hull_vals[actual_pos & model_neg]
        false_pos = e_above_hull_vals[actual_neg & model_pos]
        true_neg = e_above_hull_vals[actual_neg & model_neg]
        xlabel = r"$\Delta E_{Hull-MP}$ (eV / atom)"

    # --- histogram by model-predicted distance to convex hull
    if energy_type == "pred":
        true_pos = residuals[actual_pos & model_pos]
        false_neg = residuals[actual_pos & model_neg]
        false_pos = residuals[actual_neg & model_pos]
        true_neg = residuals[actual_neg & model_neg]
        xlabel = r"$\Delta E_{Hull-Pred}$ (eV / atom)"

    ax.hist(
        [true_pos, false_neg, false_pos, true_neg],
        bins=200,
        range=x_lim,
        alpha=0.5,
        color=["tab:green", "tab:orange", "tab:red", "tab:blue"],
        label=[
            "True Positives",
            "False Negatives",
            "False Positives",
            "True Negatives",
        ],
        stacked=True,
    )

    n_true_pos, n_false_pos, n_true_neg, n_false_neg = (
        len(true_pos),
        len(false_pos),
        len(true_neg),
        len(false_neg),
    )
    # null = (tp + fn) / (tp + tn + fp + fn)
    precision = n_true_pos / (n_true_pos + n_false_pos)

    assert n_true_pos + n_false_pos + n_true_neg + n_false_neg == len(df)

    # recall = n_true_pos / n_total_pos
    # f"Prevalence = {null:.2f}\n{precision = :.2f}\n{recall = :.2f}",
    text = f"Enrichment\nFactor = {precision/null:.3}"
    if show_mae:
        MAE = error.abs().mean()
        text += f"\n{MAE = :.3}"

    ax.text(
        0.98,
        0.98,
        text,
        fontsize=18,
        verticalalignment="top",
        horizontalalignment="right",
        transform=ax.transAxes,
    )

    ax.set(xlabel=xlabel, ylabel="Number of compounds")

    return ax


def rolling_mae_vs_hull_dist(
    df: pd.DataFrame,
    e_above_hull_col: str,
    residual_col: str = "residual",
    half_window: float = 0.02,
    bin_width: float = 0.002,
    x_lim: tuple[float, float] = (-0.2, 0.3),
    ax: plt.Axes = None,
    **kwargs: Any,
) -> plt.Axes:
    """Rolling mean absolute error as the energy to the convex hull is varied. A scale
    bar is shown for the windowing period of 40 meV per atom used when calculating
    the rolling MAE. The standard error in the mean is shaded
    around each curve. The highlighted V-shaped region shows the area in which the
    average absolute error is greater than the energy to the known convex hull. This is
    where models are most at risk of misclassifying structures.
    """
    if ax is None:
        ax = plt.gca()

    for col in (residual_col, e_above_hull_col):
        n_nans = df[col].isna().sum()
        assert n_nans == 0, f"{n_nans} NaNs in {col} column"

    is_fresh_ax = len(ax.lines) == 0

    bins = np.arange(*x_lim, bin_width)

    rolling_maes = np.zeros_like(bins)
    rolling_stds = np.zeros_like(bins)
    df = df.sort_values(by=e_above_hull_col)
    for idx, bin_center in enumerate(bins):
        low = bin_center - half_window
        high = bin_center + half_window

        mask = (df[e_above_hull_col] <= high) & (df[e_above_hull_col] > low)
        rolling_maes[idx] = df[residual_col].loc[mask].abs().mean()
        rolling_stds[idx] = scipy.stats.sem(df[residual_col].loc[mask].abs())

    ax.plot(bins, rolling_maes, **kwargs)

    ax.fill_between(
        bins, rolling_maes + rolling_stds, rolling_maes - rolling_stds, alpha=0.3
    )

    if not is_fresh_ax:
        # return earlier if all plot objects besides the line were already drawn by a
        # previous call
        return ax

    scale_bar = AnchoredSizeBar(
        ax.transData,
        2 * half_window,
        "40 meV",
        "lower left",
        pad=0.5,
        frameon=False,
        size_vertical=0.002,
    )

    ax.add_artist(scale_bar)

    ax.plot((0.05, 0.5), (0.05, 0.5), color="grey", linestyle="--", alpha=0.3)
    ax.plot((-0.5, -0.05), (0.5, 0.05), color="grey", linestyle="--", alpha=0.3)
    ax.plot((-0.05, 0.05), (0.05, 0.05), color="grey", linestyle="--", alpha=0.3)
    ax.plot((-0.1, 0.1), (0.1, 0.1), color="grey", linestyle="--", alpha=0.3)

    ax.fill_between(
        (-0.5, -0.05, 0.05, 0.5),
        (0.5, 0.5, 0.5, 0.5),
        (0.5, 0.05, 0.05, 0.5),
        color="tab:red",
        alpha=0.2,
    )

    ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
    ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

    ax.fill_between(
        (-0.05, 0, 0.05),
        (0.05, 0.05, 0.05),
        (0.05, 0, 0.05),
        color="tab:orange",
        alpha=0.2,
    )

    arrowprops = dict(facecolor="black", width=0.5, headwidth=5, headlength=5)
    ax.annotate(
        xy=(0.055, 0.05),
        xytext=(0.12, 0.05),
        arrowprops=arrowprops,
        text="Corrected\nGGA DFT\nAccuracy",
        verticalalignment="center",
        horizontalalignment="left",
    )
    ax.annotate(
        xy=(0.105, 0.1),
        xytext=(0.16, 0.1),
        arrowprops=arrowprops,
        text="GGA DFT\nAccuracy",
        verticalalignment="center",
        horizontalalignment="left",
    )

    ax.text(0, 0.13, r"$|\Delta E_{Hull-MP}| > $MAE", horizontalalignment="center")

    ax.set(xlabel=r"$\Delta E_{Hull-MP}$ (eV / atom)", ylabel="MAE (eV / atom)")

    ax.set(xlim=x_lim, ylim=(0.0, 0.14))

    return ax


def precision_recall_vs_calc_count(
    df: pd.DataFrame,
    residual_col: str = "residual",
    e_above_hull_col: str = "e_above_hull",
    stability_crit: StabilityCriterion = "energy",
    stability_threshold: float = 0,  # set stability threshold as distance to convex
    # hull in eV / atom, usually 0 or 0.1 eV
    ax: plt.Axes = None,
    label: str = None,
    intersect_lines: str | Sequence[str] = (),
    **kwargs: Any,
) -> plt.Axes:
    """Precision and recall as a function of the number of calculations performed.

    Args:
        df (pd.DataFrame): Model predictions and target energy values.
        residual_col (str, optional): Column name with residuals of model predictions,
            i.e. residual = pred - target. Defaults to "residual".
        e_above_hull_col (str, optional): Column name with convex hull distance values.
            Defaults to "e_above_hull".
        stability_crit ('energy' | 'energy+std' | 'energy-std', optional): Whether to
            use energy+/-std as stability stability_crit where std is the model
            predicted uncertainty for the energy it stipulated. Defaults to "energy".
        stability_threshold (float, optional): Max distance from convex hull before
            material is considered unstable. Defaults to 0.
        label (str, optional): Model name used to identify its liens in the legend.
            Defaults to None.
        intersect_lines (Sequence[str], optional): precision_{x,y,xy} and/or
            recall_{x,y,xy}. Defaults to (), i.e. no intersect lines.

    Returns:
        plt.Axes: The matplotlib axes object.
    """
    if ax is None:
        ax = plt.gca()

    for col in (residual_col, e_above_hull_col):
        n_nans = df[col].isna().sum()
        assert n_nans == 0, f"{n_nans} NaNs in {col} column"

    is_fresh_ax = len(ax.lines) == 0

    df = df.sort_values(by="residual")
    residuals = df[residual_col]

    if stability_crit not in get_args(StabilityCriterion):
        raise ValueError(
            f"Invalid {stability_crit=} must be one of {get_args(StabilityCriterion)}"
        )
    if "std" in stability_crit:
        # TODO column names to compute standard deviation from are currently hardcoded
        # needs to be updated when adding non-aviary models with uncertainty estimation
        var_aleatoric = (df.filter(like="_ale_") ** 2).mean(axis=1)
        var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
        std_total = (var_epistemic + var_aleatoric) ** 0.5

        if stability_crit == "energy+std":
            residuals += std_total
        elif stability_crit == "energy-std":
            residuals -= std_total

    # stability_threshold = 0.02
    stability_threshold = 0
    # stability_threshold = 0.10

    true_pos_mask = (df[e_above_hull_col] <= stability_threshold) & (
        df.residual <= stability_threshold
    )
    false_neg_mask = (df[e_above_hull_col] <= stability_threshold) & (
        df.residual > stability_threshold
    )
    false_pos_mask = (df[e_above_hull_col] > stability_threshold) & (
        df.residual <= stability_threshold
    )

    true_pos_cumsum = true_pos_mask.cumsum()

    ppv = true_pos_cumsum / (true_pos_cumsum + false_pos_mask.cumsum()) * 100
    n_true_pos = sum(true_pos_mask)
    n_false_neg = sum(false_neg_mask)
    n_total_pos = n_true_pos + n_false_neg
    tpr = true_pos_cumsum / n_total_pos * 100

    end = int(np.argmax(tpr))

    xs = np.arange(end)

    precision_curve = scipy.interpolate.interp1d(xs, ppv[:end], kind="cubic")
    rolling_recall_curve = scipy.interpolate.interp1d(xs, tpr[:end], kind="cubic")

    line_kwargs = dict(
        linewidth=4,
        markevery=[-1],
        marker="x",
        markersize=14,
        markeredgewidth=2.5,
        **kwargs,
    )
    ax.plot(xs, precision_curve(xs), linestyle="-", **line_kwargs)
    ax.plot(xs, rolling_recall_curve(xs), linestyle=":", **line_kwargs)
    ax.plot((0, 0), (0, 0), label=label, **line_kwargs)

    if intersect_lines == "all":
        intersect_lines = ("precision_xy", "recall_xy")
    if isinstance(intersect_lines, str):
        intersect_lines = [intersect_lines]
    for line_name in intersect_lines:
        try:
            line_name_map = dict(precision=precision_curve, recall=rolling_recall_curve)
            y_func = line_name_map[line_name.split("_")[0]]
        except KeyError:
            raise ValueError(
                f"Invalid {intersect_lines=}, must be one of {list(line_name_map)}"
            )
        intersect_kwargs = dict(
            linestyle=":", alpha=0.4, color=kwargs.get("color", "gray")
        )
        # Add some visual guidelines
        if "x" in line_name:
            ax.plot((0, xs[-1]), (y_func(xs[-1]), y_func(xs[-1])), **intersect_kwargs)
        if "y" in line_name:
            ax.plot((xs[-1], xs[-1]), (0, y_func(xs[-1])), **intersect_kwargs)

    if not is_fresh_ax:
        # return earlier if all plot objects besides the line were already drawn by a
        # previous call
        return ax

    ax.set(xlabel="Number of Calculations", ylabel="Precision and Recall (%)")

    ax.set(ylim=(0, 100))

    [precision] = ax.plot((0, 0), (0, 0), "black", linestyle="-")
    [recall] = ax.plot((0, 0), (0, 0), "black", linestyle=":")
    legend = ax.legend(
        [precision, recall],
        ("Precision", "Recall"),
        frameon=False,
        loc="upper right",
    )
    ax.add_artist(legend)

    return ax
