from __future__ import annotations

from typing import Any, Literal, get_args

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import scipy.interpolate
import scipy.stats
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"

StabilityCriterion = Literal["energy", "energy+std", "energy-std"]
WhichEnergy = Literal["true", "pred"]
AxLine = Literal["x", "y", "xy", ""]


# --- start global plot settings
quantity_labels = dict(
    n_atoms="Atom Count",
    n_elems="Element Count",
    crystal_sys="Crystal system",
    spg_num="Space group",
    n_wyckoff="Number of Wyckoff positions",
    n_sites="Lattice site count",
    energy_per_atom="Energy (eV/atom)",
    e_form="Formation energy (eV/atom)",
    e_above_hull="Energy above convex hull (eV/atom)",
    e_above_hull_pred="Predicted energy above convex hull (eV/atom)",
    e_above_hull_mp="Energy above MP convex hull (eV/atom)",
    e_above_hull_error="Error in energy above convex hull (eV/atom)",
    vol_diff="Volume difference (A^3)",
)
model_labels = dict(
    wren="Wren",
    wrenformer="Wrenformer",
    m3gnet="M3GNet",
    bowsr_megnet="BOWSR + MEGNet",
    cgcnn="CGCNN",
    voronoi="Voronoi",
    wbm="WBM",
    dft="DFT",
)
px.defaults.labels = quantity_labels | model_labels

pio.templates.default = "plotly_white"

# https://github.com/plotly/Kaleido/issues/122#issuecomment-994906924
# when seeing MathJax "loading" message in exported PDFs, try:
# pio.kaleido.scope.mathjax = None


plt.rc("font", size=14)
plt.rc("legend", fontsize=16, title_fontsize=16)
plt.rc("axes", titlesize=16, labelsize=16)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rc("figure", dpi=200, titlesize=16)
plt.rcParams["figure.constrained_layout.use"] = True
# --- end global plot settings


def hist_classified_stable_as_func_of_hull_dist(
    e_above_hull_pred: pd.Series,
    e_above_hull_true: pd.Series,
    std_pred: pd.Series = None,
    ax: plt.Axes = None,
    which_energy: WhichEnergy = "true",
    stability_crit: StabilityCriterion = "energy",
    stability_threshold: float = 0,
    show_threshold: bool = True,
    x_lim: tuple[float | None, float | None] = (-0.4, 0.4),
    rolling_accuracy: float = 0.02,
) -> tuple[plt.Axes, dict[str, float]]:
    """
    Histogram of the energy difference (either according to DFT ground truth [default]
    or model predicted energy) to the convex hull for materials in the WBM data set. The
    histogram is broken down into true positives, false negatives, false positives, and
    true negatives based on whether the model predicts candidates to be below the known
    convex hull. Ideally, in discovery setting a model should exhibit high recall, i.e.
    the majority of materials below the convex hull being correctly identified by the
    model.

    See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.

    Args:
        e_above_hull_pred (pd.Series): energy difference to convex hull predicted by
            model, i.e. difference between the model's predicted and true formation
            energy.
        e_above_hull_true (pd.Series): energy diff to convex hull according to DFT
            ground truth.
        std_pred (pd.Series, optional): standard deviation of the model's predicted
            formation energy.
        ax (plt.Axes, optional): matplotlib axes to plot on.
        which_energy (WhichEnergy, optional): Whether to use the true formation energy
            or the model's predicted formation energy for the histogram.
        stability_crit (StabilityCriterion, optional): Whether to add/subtract the
            model's predicted uncertainty from its energy prediction when measuring
            predicted stability.
        stability_threshold (float, optional): set stability threshold as distance to
            convex hull in eV/atom, usually 0 or 0.1 eV.
        show_threshold (bool, optional): Whether to plot stability threshold as dashed
            vertical line.
        x_lim (tuple[float | None, float | None]): x-axis limits.
        rolling_accuracy (float): Rolling accuracy window size in eV / atom. Set to 0 to
            disable. Defaults to 0.01.

    Returns:
        tuple[plt.Axes, dict[str, float]]: plot axes and classification metrics

    NOTE this figure plots hist bars separately which causes aliasing in pdf. Can be
    fixed in Inkscape or similar by merging regions by color.
    """
    ax = ax or plt.gca()

    if stability_crit not in get_args(StabilityCriterion):
        raise ValueError(
            f"Invalid {stability_crit=} must be one of {get_args(StabilityCriterion)}"
        )

    test = e_above_hull_pred + e_above_hull_true
    if stability_crit == "energy+std":
        test += std_pred
    elif stability_crit == "energy-std":
        test -= std_pred

    # --- histogram of DFT-computed distance to convex hull
    if which_energy == "true":
        actual_pos = e_above_hull_true <= stability_threshold
        actual_neg = e_above_hull_true > stability_threshold
        model_pos = test <= stability_threshold
        model_neg = test > stability_threshold

        n_true_pos = len(e_above_hull_true[actual_pos & model_pos])
        n_false_neg = len(e_above_hull_true[actual_pos & model_neg])

        n_total_pos = n_true_pos + n_false_neg
        null = n_total_pos / len(e_above_hull_true)

        true_pos = e_above_hull_true[actual_pos & model_pos]
        false_neg = e_above_hull_true[actual_pos & model_neg]
        false_pos = e_above_hull_true[actual_neg & model_pos]
        true_neg = e_above_hull_true[actual_neg & model_neg]
        xlabel = r"$\Delta E_{Hull-MP}$ (eV / atom)"

    # --- histogram of model-predicted distance to convex hull
    if which_energy == "pred":
        true_pos = e_above_hull_pred[actual_pos & model_pos]
        false_neg = e_above_hull_pred[actual_pos & model_neg]
        false_pos = e_above_hull_pred[actual_neg & model_pos]
        true_neg = e_above_hull_pred[actual_neg & model_neg]
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

    n_true_pos, n_false_pos, n_true_neg, n_false_neg = map(
        len, (true_pos, false_pos, true_neg, false_neg)
    )
    # null = (tp + fn) / (tp + tn + fp + fn)
    precision = n_true_pos / (n_true_pos + n_false_pos)

    # assert (n_all := n_true_pos + n_false_pos + n_true_neg + n_false_neg) == len(
    #     e_above_hull_true
    # ), f"{n_all} != {len(e_above_hull_true)}"

    ax.set(xlabel=xlabel, ylabel="Number of compounds", xlim=x_lim)

    if rolling_accuracy:
        # add moving average of the accuracy (computed within 20 meV/atom intervals) as
        # a function of ΔHd,MP is shown as a blue line (right axis)
        ax_acc = ax.twinx()
        ax_acc.set_ylabel("Accuracy", color="darkblue")
        ax_acc.tick_params(labelcolor="darkblue")
        ax_acc.set(ylim=(0, 1))

        # --- moving average of the accuracy
        # compute accuracy within 20 meV/atom intervals
        bins = np.arange(x_lim[0], x_lim[1], rolling_accuracy)
        bin_counts = np.histogram(e_above_hull_true, bins)[0]
        bin_true_pos = np.histogram(true_pos, bins)[0]
        bin_true_neg = np.histogram(true_neg, bins)[0]

        # compute accuracy
        bin_accuracies = (bin_true_pos + bin_true_neg) / bin_counts
        # plot accuracy
        ax_acc.plot(
            bins[:-1],
            bin_accuracies,
            color="tab:blue",
            label="Accuracy",
            linewidth=3,
        )
        # ax2.fill_between(
        #     bin_centers,
        #     bin_accuracy - bin_accuracy_std,
        #     bin_accuracy + bin_accuracy_std,
        #     color="tab:blue",
        #     alpha=0.2,
        # )

    if show_threshold:
        ax.axvline(
            stability_threshold,
            color="k",
            linestyle="--",
            label="Stability Threshold",
        )

    recall = n_true_pos / n_total_pos

    return ax, {
        "enrichment": precision / null,
        "precision": precision,
        "recall": recall,
        "prevalence": null,
        "accuracy": (n_true_pos + n_true_neg)
        / (n_true_pos + n_true_neg + n_false_pos + n_false_neg),
        "f1": 2 * (precision * recall) / (precision + recall),
    }


def rolling_mae_vs_hull_dist(
    e_above_hull_true: pd.Series,
    e_above_hull_pred: pd.Series,
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
    ax = ax or plt.gca()

    for series in (e_above_hull_pred, e_above_hull_true):
        n_nans = series.isna().sum()
        assert n_nans == 0, f"{n_nans:,} NaNs in {series.name}"

    is_fresh_ax = len(ax.lines) == 0

    bins = np.arange(*x_lim, bin_width)

    rolling_maes = np.zeros_like(bins)
    rolling_stds = np.zeros_like(bins)
    for idx, bin_center in enumerate(bins):
        low = bin_center - half_window
        high = bin_center + half_window

        mask = (e_above_hull_true <= high) & (e_above_hull_true > low)
        rolling_maes[idx] = e_above_hull_pred.loc[mask].abs().mean()
        rolling_stds[idx] = scipy.stats.sem(e_above_hull_pred.loc[mask].abs())

    kwargs = dict(linewidth=3) | kwargs
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

    # shrink=0.05 means cut off 5% length from both sides of arrow line
    arrowprops = dict(
        facecolor="black", width=0.5, headwidth=5, headlength=5, shrink=0.05
    )
    ax.annotate(
        xy=(0.05, 0.05),
        xytext=(0.2, 0.05),
        text="Corrected\nGGA DFT\nAccuracy",
        arrowprops=arrowprops,
        verticalalignment="center",
    )
    ax.annotate(
        xy=(0.1, 0.1),
        xytext=(0.2, 0.1),
        text="GGA DFT\nAccuracy",
        arrowprops=arrowprops,
        verticalalignment="center",
    )

    ax.text(0, 0.13, r"$|\Delta E_{Hull-MP}| > $MAE", horizontalalignment="center")
    ax.set(xlabel=r"$\Delta E_{Hull-MP}$ (eV / atom)", ylabel="MAE (eV / atom)")
    ax.set(xlim=x_lim, ylim=(0.0, 0.14))

    return ax


def cumulative_clf_metric(
    e_above_hull_error: pd.Series,
    e_above_hull_true: pd.Series,
    metric: Literal["precision", "recall"],
    std_pred: pd.Series = None,
    stability_crit: StabilityCriterion = "energy",
    stability_threshold: float = 0,  # set stability threshold as distance to convex
    # hull in eV / atom, usually 0 or 0.1 eV
    ax: plt.Axes = None,
    label: str = None,
    project_end_point: AxLine = "xy",
    show_optimal: bool = False,
    **kwargs: Any,
) -> plt.Axes:
    """Precision and recall as a function of the number of included materials sorted
    by model-predicted distance to the convex hull, i.e. materials predicted most stable
    enter the precision and recall calculation first. The curves end when all materials
    predicted stable are included.

    Args:
        df (pd.DataFrame): Model predictions and target energy values.
        e_above_hull_error (str, optional): Column name with residuals of model
            predictions, i.e. residual = pred - target. Defaults to "residual".
        e_above_hull_true (str, optional): Column name with convex hull distance values.
            Defaults to "e_above_hull".
        metric ('precision' | 'recall', optional): Metric to plot.
        stability_crit ('energy' | 'energy+std' | 'energy-std', optional): Whether to
            use energy+/-std as stability stability_crit where std is the model
            predicted uncertainty for the energy it stipulated. Defaults to "energy".
        stability_threshold (float, optional): Max distance from convex hull before
            material is considered unstable. Defaults to 0.
        label (str, optional): Model name used to identify its liens in the legend.
            Defaults to None.
        project_end_point ('x' | 'y' | 'xy' | '', optional): Defaults to '', i.e. no
            axis projection lines.
        show_optimal (bool, optional): Whether to plot the optimal precision/recall
            line. Defaults to False.

    Returns:
        plt.Axes: The matplotlib axes object.
    """
    ax = ax or plt.gca()

    for series in (e_above_hull_error, e_above_hull_true):
        n_nans = series.isna().sum()
        assert n_nans == 0, f"{n_nans:,} NaNs in {series.name}"

    e_above_hull_error = e_above_hull_error.sort_values()
    e_above_hull_true = e_above_hull_true.loc[e_above_hull_error.index]

    if stability_crit not in get_args(StabilityCriterion):
        raise ValueError(
            f"Invalid {stability_crit=} must be one of {get_args(StabilityCriterion)}"
        )
    if stability_crit == "energy+std":
        e_above_hull_error += std_pred
    elif stability_crit == "energy-std":
        e_above_hull_error -= std_pred

    true_pos_mask = (e_above_hull_true <= stability_threshold) & (
        e_above_hull_error <= stability_threshold
    )
    false_neg_mask = (e_above_hull_true <= stability_threshold) & (
        e_above_hull_error > stability_threshold
    )
    false_pos_mask = (e_above_hull_true > stability_threshold) & (
        e_above_hull_error <= stability_threshold
    )

    true_pos_cumsum = true_pos_mask.cumsum()

    # precision aka positive predictive value (PPV)
    precision = true_pos_cumsum / (true_pos_cumsum + false_pos_mask.cumsum()) * 100
    n_true_pos = sum(true_pos_mask)
    n_false_neg = sum(false_neg_mask)
    n_total_pos = n_true_pos + n_false_neg
    true_pos_rate = true_pos_cumsum / n_total_pos * 100

    end = int(np.argmax(true_pos_rate))
    xs = np.arange(end)

    ys_raw = dict(precision=precision, recall=true_pos_rate)[metric]
    y_interp = scipy.interpolate.interp1d(xs, ys_raw[:end], kind="cubic")
    ys = y_interp(xs)

    line_kwargs = dict(
        linewidth=2, markevery=[-1], marker="x", markersize=14, markeredgewidth=2.5
    )
    ax.plot(xs, ys, **line_kwargs | kwargs)
    ax.text(
        xs[-1],
        ys[-1],
        label,
        color=kwargs.get("color"),
        verticalalignment="bottom",
        rotation=30,
        bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"),
    )

    # add some visual guidelines
    intersect_kwargs = dict(linestyle=":", alpha=0.4, color=kwargs.get("color"))
    if "x" in project_end_point:
        ax.plot((0, xs[-1]), (ys[-1], ys[-1]), **intersect_kwargs)
    if "y" in project_end_point:
        ax.plot((xs[-1], xs[-1]), (0, ys[-1]), **intersect_kwargs)

    ax.set(ylim=(0, 100), ylabel=f"{metric.title()} (%)")

    # optimal recall line finds all stable materials without any false positives
    # can be included to confirm all models start out of with near optimal recall
    # and to see how much each model overshoots total n_stable
    n_below_hull = sum(e_above_hull_true < 0)
    if show_optimal:
        ax.plot(
            [0, n_below_hull],
            [0, 100],
            color="green",
            linestyle="dashed",
            linewidth=1,
            label=f"Optimal {metric.title()}",
        )
        ax.text(
            n_below_hull,
            100,
            label,
            color=kwargs.get("color"),
            verticalalignment="top",
            rotation=-30,
            bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"),
        )

    return ax
