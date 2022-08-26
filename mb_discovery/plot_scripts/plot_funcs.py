from __future__ import annotations

from typing import Literal, Sequence

import matplotlib.pyplot as plt
import pandas as pd


__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"


plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
plt.rc("font", size=16)


def hist_classified_stable_as_func_of_hull_dist(
    df: pd.DataFrame,
    target_col: str,
    pred_cols: Sequence[str],
    e_above_hull_col: str,
    ax: plt.Axes = None,
    energy_type: Literal["true", "pred"] = "true",
    criterion: Literal["energy", "std", "neg_std"] = "energy",
    show_mae: bool = False,
    stability_thresh: float = 0,  # set stability threshold as distance to convex hull
    # in eV / atom, usually 0 or 0.1 eV
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
    mean = error + e_above_hull_vals

    if criterion == "energy":
        test = mean
    elif "std" in criterion:
        # TODO column names to compute standard deviation from are currently hardcoded
        # needs to be updated when adding non-aviary models with uncertainty estimation
        var_aleatoric = (df.filter(like="_ale_") ** 2).mean(axis=1)
        var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
        std_total = (var_epistemic + var_aleatoric) ** 0.5

        if criterion == "std":
            test += std_total
        elif criterion == "neg_std":
            test -= std_total

    # --- histogram by DFT-computed distance to convex hull
    if energy_type == "true":
        actual_pos = e_above_hull_vals <= stability_thresh
        actual_neg = e_above_hull_vals > stability_thresh
        model_pos = test <= stability_thresh
        model_neg = test > stability_thresh

        n_true_pos = len(e_above_hull_vals[actual_pos & model_pos])
        n_false_neg = len(e_above_hull_vals[actual_pos & model_neg])

        n_total_pos = n_true_pos + n_false_neg
        null = n_total_pos / len(e_above_hull_vals)

        true_pos = e_above_hull_vals[actual_pos & model_pos]
        false_neg = e_above_hull_vals[actual_pos & model_neg]
        false_pos = e_above_hull_vals[actual_neg & model_pos]
        true_neg = e_above_hull_vals[actual_neg & model_neg]
        xlabel = r"$\Delta E_{Hull-MP}$ / eV per atom"

    # --- histogram by model-predicted distance to convex hull
    if energy_type == "pred":
        true_pos = mean[actual_pos & model_pos]
        false_neg = mean[actual_pos & model_neg]
        false_pos = mean[actual_neg & model_pos]
        true_neg = mean[actual_neg & model_neg]
        xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

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
