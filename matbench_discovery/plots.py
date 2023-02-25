from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
import scipy.interpolate
import scipy.stats
import wandb
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from tqdm import tqdm

from matbench_discovery.metrics import classify_stable

__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"

Backend = Literal["matplotlib", "plotly"]

# --- start global plot settings
quantity_labels = dict(
    n_atoms="Atom Count",
    n_elems="Element Count",
    crystal_sys="Crystal system",
    spg_num="Space group",
    n_wyckoff="Number of Wyckoff positions",
    n_sites="Lattice site count",
    energy_per_atom="Energy (eV/atom)",
    e_form="Actual E<sub>form</sub> (eV/atom)",
    e_above_hull="E<sub>above hull</sub> (eV/atom)",
    e_above_hull_mp2020_corrected_ppd_mp="Actual E<sub>above hull</sub> (eV/atom)",
    e_above_hull_pred="Predicted E<sub>above hull</sub> (eV/atom)",
    e_above_hull_mp="E<sub>above MP hull</sub> (eV/atom)",
    e_above_hull_error="Error in E<sub>above hull</sub> (eV/atom)",
    vol_diff="Volume difference (A^3)",
    e_form_per_atom_mp2020_corrected="Actual E<sub>form</sub> (eV/atom)",
    e_form_per_atom_pred="Predicted E<sub>form</sub> (eV/atom)",
    material_id="Material ID",
    band_gap="Band gap (eV)",
    formula="Formula",
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

# color list https://plotly.com/python-api-reference/generated/plotly.graph_objects.layout
colorway = (
    "lightseagreen",
    "orange",
    "lightsalmon",
    "dodgerblue",
)
clf_labels = ("True Positive", "False Negative", "False Positive", "True Negative")
clf_colors = ("lightseagreen", "orange", "lightsalmon", "dodgerblue")
clf_color_map = dict(zip(clf_labels, clf_colors))

global_layout = dict(
    # colorway=px.colors.qualitative.Pastel,
    # colorway=colorway,
    margin=dict(l=30, r=20, t=60, b=20),
    paper_bgcolor="rgba(0,0,0,0)",
    # plot_bgcolor="rgba(0,0,0,0)",
    font_size=13,
    # increase legend marker size and make background transparent
    legend=dict(itemsizing="constant", bgcolor="rgba(0, 0, 0, 0)"),
)
pio.templates["global"] = dict(layout=global_layout)
pio.templates.default = "plotly_dark+global"

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


def hist_classified_stable_vs_hull_dist(
    df: pd.DataFrame,
    each_true_col: str,
    each_pred_col: str,
    ax: plt.Axes = None,
    which_energy: Literal["true", "pred"] = "true",
    stability_threshold: float | None = 0,
    x_lim: tuple[float | None, float | None] = (-0.7, 0.7),
    rolling_acc: float | None = 0.02,
    backend: Backend = "plotly",
    y_label: str = "Number of materials",
    clf_labels: Sequence[str] = (
        "True Positive",
        "False Negative",
        "False Positive",
        "True Negative",
    ),
    **kwargs: Any,
) -> plt.Axes | go.Figure:
    """Histogram of the energy difference (either according to DFT ground truth [default]
    or model predicted energy) to the convex hull for materials in the WBM data set. The
    histogram is broken down into true positives, false negatives, false positives, and
    true negatives based on whether the model predicts candidates to be below the known
    convex hull. Ideally, in discovery setting a model should exhibit high recall, i.e.
    the majority of materials below the convex hull being correctly identified by the
    model.

    See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.

    Args:
        df (pd.DataFrame): Data frame containing true and predicted hull distances.
        each_true_col (str): Name of column with energy above convex hull according to DFT
            ground truth (in eV / atom).
        each_pred_col (str): Name of column with energy above convex hull predicted by model
            (in eV / atom). Same as true energy to convex hull plus predicted minus true
            formation energy.
        ax (plt.Axes, optional): matplotlib axes to plot on.
        which_energy ('true' | 'pred', optional): Whether to use the true (DFT) hull
            distance or the model's predicted hull distance for the histogram.
        stability_threshold (float, optional): set stability threshold as distance to
            convex hull in eV/atom, usually 0 or 0.1 eV.
        x_lim (tuple[float | None, float | None]): x-axis limits.
        rolling_acc (float): Rolling accuracy window size in eV / atom. Set to None
            or 0 to disable. Defaults to 0.02, meaning 20 meV / atom.
        backend ('matplotlib' | 'plotly'], optional): Which plotting engine to use.
            Changes the return type. Defaults to 'plotly'.
        y_label (str, optional): y-axis label. Defaults to "Number of materials".
        clf_labels (list[str], optional): Labels for the four classification categories.
            Defaults to ["True Positive", "False Negative", "False Positive", "True
            Negative"].
        kwargs: Additional keyword arguments passed to the ax.hist() or px.histogram()
            depending on backend.

    Returns:
        tuple[plt.Axes, dict[str, float]]: plot axes and classification metrics

    NOTE this figure plots hist bars separately which causes aliasing in pdf. Can be
    fixed in Inkscape or similar by merging regions by color.
    """
    each_true = df[each_true_col]
    each_pred = df[each_pred_col]
    x_col = each_true_col if which_energy == "true" else each_pred_col
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        each_true, each_pred, stability_threshold
    )

    # toggle between histogram of DFT-computed or model-predicted distance to convex hull
    e_above_hull = df[x_col]
    eah_true_pos = e_above_hull[true_pos]
    eah_true_neg = e_above_hull[true_neg]
    # eah_false_neg = e_above_hull[false_neg]
    # eah_false_pos = e_above_hull[false_pos]
    # n_true_pos, n_false_pos, n_true_neg, n_false_neg = map(
    #     sum, (true_pos, false_pos, true_neg, false_neg)
    # )

    clf = np.array(clf_labels)[
        true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
    ]

    df[(clf_col := "classified")] = clf

    kwds: dict[str, Any] = (
        dict(
            barmode="stack",
            range_x=x_lim,
            color_discrete_map=clf_color_map,
            x=x_col,
            color=clf_col,
        )
        if backend == "plotly"
        else dict(
            column=[x_col],
            legend=False,
            ax=ax,
            xlim=x_lim,
            # color=df[clf_col],
        )
    )

    fig = df.plot.hist(backend=backend, **(kwds | kwargs))

    if backend == "matplotlib":
        # ax.hist(
        #     [eah_true_pos, eah_false_neg, eah_false_pos, eah_true_neg],
        #     bins=200,
        #     range=x_lim,
        #     alpha=0.5,
        #     color=["tab:green", "tab:orange", "tab:red", "tab:blue"],
        #     label=clf_labels,
        #     stacked=True,
        #     **kwargs,
        # )
        xlabel = dict(
            true=r"$E_\mathrm{above\ hull}\;\mathrm{(eV / atom)}$",
            pred=r"$E_\mathrm{above\ hull\ pred}\;\mathrm{(eV / atom)}$",
        )[which_energy]

        if stability_threshold is not None:
            for ax in [fig] if isinstance(fig, plt.Axes) else fig.flat:
                ax.set(xlabel=xlabel, ylabel=y_label, xlim=x_lim)
                label = "Stability Threshold"
                ax.axvline(
                    stability_threshold, color="black", linestyle="--", label=label
                )

    if backend == "plotly":
        fig.update_layout(legend=dict(title=None, y=0.5, xanchor="right", x=1))

        if stability_threshold is not None:
            anno = dict(text="Stability threshold", xanchor="center", yanchor="bottom")
            fig.add_vline(stability_threshold, line=dict(dash="dash"), annotation=anno)

    if rolling_acc:
        # add moving average of the accuracy computed within given window
        # as a function of e_above_hull shown as blue line (right axis)

        # --- moving average of the accuracy
        # compute rolling accuracy in rolling_acc-sized intervals
        bins = np.arange(each_true.min(), each_true.max(), rolling_acc)
        bin_counts = np.histogram(each_true, bins)[0]
        bin_true_pos = np.histogram(eah_true_pos, bins)[0]
        bin_true_neg = np.histogram(eah_true_neg, bins)[0]

        # compute accuracy (handling division by zero)
        bin_accuracies = np.divide(
            bin_true_pos + bin_true_neg,
            bin_counts,
            out=np.zeros_like(bin_counts, dtype=float),
            where=bin_counts != 0,
        )

        if backend == "matplotlib":
            for ax in fig.flat if isinstance(fig, np.ndarray) else [fig]:
                ax_acc = ax.twinx()
                ax_acc.set_ylabel("Rolling Accuracy", color="darkblue")
                ax_acc.tick_params(labelcolor="darkblue")
                ax_acc.set(ylim=(0, 1.1))
                # plot accuracy
                ax_acc.plot(
                    bins[:-1],
                    bin_accuracies,
                    color="tab:blue",
                    label="Accuracy",
                    linewidth=3,
                )
                # ax_acc.fill_between(
                #     bin_centers,
                #     bin_accuracies - bin_accuracy_std,
                #     bin_accuracies + bin_accuracy_std,
                #     color="tab:blue",
                #     alpha=0.2,
                # )
        else:
            style = dict(color="orange")
            title = "Rolling Accuracy"
            # add accuracy line on a separate y-axis on the right side of the plot
            fig.update_layout(
                yaxis2=dict(overlaying="y", side="right", range=[0, 1], tickfont=style),
                # title = dict(text=title, font=style)
            )
            fig.add_scatter(
                x=bins, y=bin_accuracies, name=title, line=style, yaxis="y2"
            )

    return fig


def rolling_mae_vs_hull_dist(
    e_above_hull_true: pd.Series,
    e_above_hull_errors: pd.DataFrame | dict[str, pd.Series],
    df_rolling_err: pd.DataFrame = None,
    df_err_std: pd.DataFrame = None,
    window: float = 0.02,
    bin_width: float = 0.005,
    x_lim: tuple[float, float] = (-0.2, 0.2),
    y_lim: tuple[float, float] = (0, 0.2),
    backend: Backend = "plotly",
    y_label: str = "rolling MAE (eV/atom)",
    just_plot_lines: bool = False,
    with_sem: bool = True,
    show_dft_acc: bool = False,
    show_dummy_mae: bool = False,
    **kwargs: Any,
) -> plt.Axes | go.Figure:
    """Rolling mean absolute error as the energy to the convex hull is varied. A scale
    bar is shown for the windowing period of 40 meV per atom used when calculating the
    rolling MAE. The standard error in the mean is shaded around each curve. The
    highlighted V-shaped region shows the area in which the average absolute error is
    greater than the energy to the known convex hull. This is where models are most at
    risk of misclassifying structures.

    Args:
        e_above_hull_true (pd.Series): Distance to convex hull according to DFT
            ground truth (in eV / atom).
        e_above_hull_errors (pd.DataFrame | dict[str, pd.Series]): Error in
            model-predicted distance to convex hull, i.e. actual hull distance minus
            predicted hull distance (in eV / atom).
        df_rolling_err (pd.DataFrame, optional): Cached rolling MAE(s) as returned by
            previous call to this function. Defaults to None.
        df_err_std (pd.DataFrame, optional): Cached standard error in the mean of the
            rolling MAE(s) as returned by prev. call to this function. Defaults to None.
        window (float, optional): Rolling MAE averaging window. Defaults to 0.02 (20
            meV/atom)
        bin_width (float, optional): Density of line points (more points the
            smaller). Defaults to 0.002.
        x_lim (tuple[float, float], optional): x-axis range. Defaults to (-0.2, 0.3).
        y_lim (tuple[float, float], optional): y-axis range. Defaults to (0.0, 0.14).
        backend ('matplotlib' | 'plotly'], optional): Which plotting engine to use.
            Changes the return type. Defaults to 'plotly'.
        y_label (str, optional): y-axis label. Defaults to "rolling MAE (eV/atom)".
        just_plot_line (bool, optional): If True, plot only the rolling MAE, no shapes
            and annotations. Also won't plot the standard error in the mean. Defaults
            to False.
        just_plot_lines (bool, optional): If True, plot only the rolling MAE, no shapes
            and annotations. Also won't plot the standard error in the mean. Defaults
            to False.
        with_sem (bool, optional): If True, plot the standard error of the mean as
            shaded area around the rolling MAE. Defaults to True.
        show_dft_acc (bool, optional): If True, change color of the triangle of peril's tip
            and annotate it with 'Corrected GGA Accuracy' at rolling MAE of 25 meV/atom.
            Defaults to False.
        show_dummy_mae (bool, optional): If True, plot a line at the dummy MAE of always
            predicting the target mean.
        **kwargs: Additional keyword arguments to pass to df.plot().

    Returns:
        tuple[plt.Axes | go.Figure, pd.DataFrame, pd.DataFrame]: matplotlib Axes or
        plotly
            Figure depending on backend, followed by two dataframes containing the
            rolling error for each column in e_above_hull_errors and the rolling
            standard error in the mean.
    """
    bins = np.arange(*x_lim, bin_width)
    models = list(e_above_hull_errors)

    if df_rolling_err is None or df_err_std is None:
        df_rolling_err = pd.DataFrame(columns=models, index=bins)
        df_err_std = df_rolling_err.copy()

        for model in (pbar := tqdm(models, desc="Calculating rolling MAE")):
            pbar.set_postfix_str(model)
            for idx, bin_center in enumerate(bins):
                low = bin_center - window
                high = bin_center + window

                mask = (e_above_hull_true <= high) & (e_above_hull_true > low)

                bin_mae = e_above_hull_errors[model].loc[mask].abs().mean()
                df_rolling_err[model].iloc[idx] = bin_mae

                # drop NaNs to avoid error, scipy doesn't ignore NaNs
                each_std = scipy.stats.sem(
                    e_above_hull_errors[model].loc[mask].dropna().abs()
                )
                df_err_std[model].iloc[idx] = each_std
    else:
        print("Using pre-calculated rolling MAE")

    ax = df_rolling_err.plot(backend=backend, **kwargs)

    if just_plot_lines:
        # return earlier if all plot objects besides the line were already drawn by a
        # previous call
        return ax, df_rolling_err, df_err_std

    # DFT accuracy at 25 meV/atom for relative difference of e_above_hull for chemically
    # similar systems which is lower than formation energy error due to systematic error
    # cancellation among similar chemistries, supporting ref:
    href = "https://doi.org/10.1103/PhysRevB.85.155208"
    dft_acc = 0.025

    window_bar_anno = f"rolling window={2 * window * 1000:.0f} meV"
    dummy_mae = (e_above_hull_true - e_above_hull_true.mean()).abs().mean()
    dummy_mae_text = f"dummy MAE = {dummy_mae:.2f} eV/atom"

    if backend == "matplotlib":
        # assert df_rolling_err.isna().sum().sum() == 0, "NaNs in df_rolling_err"
        # assert df_err_std.isna().sum().sum() == 0, "NaNs in df_err_std"
        # for model in df_rolling_err if with_sem else []:
        #     ax.fill_between(
        #         bins,
        #         df_rolling_err[model] + df_err_std[model],
        #         df_rolling_err[model] - df_err_std[model],
        #         alpha=0.3,
        #     )

        scale_bar = AnchoredSizeBar(
            ax.transData,
            window,
            window_bar_anno,
            "lower left",
            pad=0.5,
            frameon=False,
            size_vertical=0.002,
        )
        # indicate size of MAE averaging window
        ax.add_artist(scale_bar)

        ax.fill_between(
            (-1, -dft_acc, dft_acc, 1) if show_dft_acc else (-1, 0, 1),
            (1, 1, 1, 1) if show_dft_acc else (1, 1, 1),
            (1, dft_acc, dft_acc, 1) if show_dft_acc else (1, 0, 1),
            color="tab:red",
            alpha=0.2,
        )

        if show_dft_acc:
            ax.fill_between(
                (-dft_acc, 0, dft_acc),
                (dft_acc, dft_acc, dft_acc),
                (dft_acc, 0, dft_acc),
                color="tab:orange",
                alpha=0.2,
            )
            # shrink=0.1 means cut off 10% length from both sides of arrow line
            arrowprops = dict(
                facecolor="black", width=0.5, headwidth=5, headlength=5, shrink=0.1
            )
            ax.annotate(
                xy=(-dft_acc, dft_acc),
                xytext=(-2 * dft_acc, dft_acc),
                text="Corrected GGA\nAccuracy",
                arrowprops=arrowprops,
                verticalalignment="center",
                horizontalalignment="right",
            )

        ax.axhline(dummy_mae, color="tab:blue", linestyle="--", linewidth=0.5)
        ax.text(dummy_mae, 0.1, dummy_mae_text)

        ax.text(
            0, 0.13, r"MAE > $|E_\mathrm{above\ hull}|$", horizontalalignment="center"
        )
        ax.set(xlabel=r"$E_\mathrm{above\ hull}$ (eV/atom)", ylabel=y_label)
        ax.set(xlim=x_lim, ylim=y_lim)

    elif backend == "plotly":
        for idx, model in enumerate(df_rolling_err if with_sem else []):
            # set legendgroup to model name so SEM shading toggles with model curve
            ax.data[idx].legendgroup = model
            # set SEM area to same color as model curve
            ax.add_scatter(
                x=list(bins) + list(bins)[::-1],  # bins, then bins reversed
                # upper, then lower reversed
                y=list(df_rolling_err[model] + 3 * df_err_std[model])
                + list(df_rolling_err[model] - 3 * df_err_std[model])[::-1],
                mode="lines",
                line=dict(color="white", width=0),
                fill="toself",
                legendgroup=model,
                fillcolor=ax.data[idx].line.color,
                opacity=0.4,
                showlegend=False,
            )

        ax.layout.legend.update(
            title="",
            x=0,
            y=0,
            # xanchor="right",
            yanchor="bottom",
            title_font=dict(size=15),
        )
        ax.layout.xaxis.title.text = "E<sub>above MP hull</sub> (eV/atom)"
        ax.layout.yaxis.title.text = "rolling MAE (eV/atom)"
        ax.update_xaxes(range=x_lim)
        ax.update_yaxes(range=y_lim)
        # exclude from hover tooltip
        scatter_kwds = dict(
            fill="toself", opacity=0.2, hoverinfo="skip", showlegend=False
        )
        peril_cone_anno = "MAE > |E<sub>above hull</sub>|"
        ax.add_scatter(
            x=(-1, -dft_acc, dft_acc, 1) if show_dft_acc else (-1, 0, 1),
            y=(1, dft_acc, dft_acc, 1) if show_dft_acc else (1, 0, 1),
            name=peril_cone_anno,
            fillcolor="red",
            **scatter_kwds,
        )
        ax.add_annotation(
            x=0,
            y=0.8,
            text=peril_cone_anno,
            showarrow=False,
            yref="paper",
        )

        if show_dummy_mae:
            ax.add_hline(
                y=dummy_mae,
                line=dict(dash="dash", width=0.5),
                annotation_text=dummy_mae_text,
            )

        if show_dft_acc:
            ax.add_scatter(
                x=(-dft_acc, dft_acc, 0, -dft_acc),
                y=(dft_acc, dft_acc, 0, dft_acc),
                name="MAE < |Corrected GGA error|",
                fillcolor="red",
                **scatter_kwds,
            )
            ax.add_annotation(
                x=-dft_acc,
                y=dft_acc,
                text=f"<a {href=}>Corrected GGA Accuracy<br>for rel. Energy</a> "
                "[<a href='#hautier_accuracy_2012' target='_self'>ref</a>]",
                showarrow=True,
                xshift=-10,
                arrowhead=2,
                ax=-4 * dft_acc,
                ay=2 * dft_acc,
                axref="x",
                ayref="y",
            )

        ax.data = ax.data[::-1]  # bring px.line() to front
        # plot rectangle to indicate MAE window size
        x0, y0 = x_lim[1] - 0.01, y_lim[0] + 0.01
        ax.add_annotation(
            x=x0 - window,
            y=y0,
            text=window_bar_anno,
            showarrow=False,
            xshift=-4,
            yshift=-6,
            yanchor="bottom",
            xanchor="right",
        )
        ax.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x0 - window,
            y1=y0 + window / 5,
        )

    return ax, df_rolling_err, df_err_std


def cumulative_precision_recall(
    e_above_hull_true: pd.Series,
    df_preds: pd.DataFrame,
    metrics: Sequence[str] = ("Cumulative Precision", "Cumulative Recall"),
    stability_threshold: float = 0,  # set stability threshold as distance to convex
    # hull in eV / atom, usually 0 or 0.1 eV
    project_end_point: Literal["x", "y", "xy", ""] = "xy",
    optimal_recall: str | None = "Optimal Recall",
    show_n_stable: bool = True,
    backend: Backend = "plotly",
    **kwargs: Any,
) -> tuple[plt.Figure | go.Figure, pd.DataFrame]:
    """Create 2 subplots side-by-side with cumulative precision and recall curves for
    all models starting with materials predicted most stable, adding the next material,
    recomputing the cumulative metrics, adding the next most stable material and so on
    until each model no longer predicts the material to be stable. Again, materials
    predicted more stable enter the precision and recall calculation sooner. Different
    models predict different number of materials to be stable. Hence the curves end at
    different points.

    Args:
        e_above_hull_true (pd.Series): Distance to convex hull according to DFT
            ground truth (in eV / atom).
        df_preds (pd.DataFrame): Distance to convex hull predicted by models, one column
            per model (in eV / atom). Same as true energy to convex hull plus predicted
            minus true formation energy.
        metrics (Sequence[str], optional): Which metrics to plot. Defaults to
            ('Precision', 'Recall'). Also accepts 'F1'.
        stability_threshold (float, optional): Max distance above convex hull before
            material is considered unstable. Defaults to 0.
        project_end_point ('x' | 'y' | 'xy' | '', optional): Whether to project end
        points of precision and recall curves to the x/y axis. Defaults to '', i.e. no
            axis projection lines.
        optimal_recall (str | None, optional): Label for the optimal recall line. Defaults
            to 'Optimal Recall'. Set to None to not plot the line.
        show_n_stable (bool, optional): Whether to show a horizontal line at the true number
            of stable materials. Defaults to True.
        backend ('matplotlib' | 'plotly'], optional): Which plotting engine to use.
            Changes the return type. Defaults to 'plotly'.
        **kwargs: Keyword arguments passed to df.plot().

    Returns:
        tuple[plt.Figure | go.Figure, pd.DataFrame]: The matplotlib/plotly figure and
            dataframe of cumulative metrics for each model.
    """
    factory = lambda: pd.DataFrame(index=range(len(e_above_hull_true)))
    dfs: dict[str, pd.DataFrame] = defaultdict(factory)

    for model_name in df_preds:
        each_pred = df_preds[model_name].sort_values()
        # sort targets by model ranking
        each_true = e_above_hull_true.loc[each_pred.index]

        true_pos_cum, false_neg_cum, false_pos_cum, _true_neg_cum = map(
            np.cumsum, classify_stable(each_true, each_pred, stability_threshold)
        )

        # precision aka positive predictive value (PPV)
        precision_cum = true_pos_cum / (true_pos_cum + false_pos_cum)
        n_total_pos = true_pos_cum[-1] + false_neg_cum[-1]
        recall_cum = true_pos_cum / n_total_pos  # aka true_pos_rate aka sensitivity
        # cumulative F1 score
        f1_cum = 2 * (precision_cum * recall_cum) / (precision_cum + recall_cum)

        end = int(np.argmax(recall_cum))
        xs = np.arange(end)
        prec_interp = scipy.interpolate.interp1d(xs, precision_cum[:end], kind="cubic")
        recall_interp = scipy.interpolate.interp1d(xs, recall_cum[:end], kind="cubic")
        f1_interp = scipy.interpolate.interp1d(xs, f1_cum[:end], kind="cubic")

        dfs["Cumulative Precision"][model_name] = pd.Series(prec_interp(xs))
        dfs["Cumulative Recall"][model_name] = pd.Series(recall_interp(xs))
        dfs["Cumulative F1"][model_name] = pd.Series(f1_interp(xs))

    for key, df in dfs.items():
        # drop all-NaN rows so plotly plot x-axis only extends to largest number of
        # predicted materials by any model
        df.dropna(how="all", inplace=True)
        # will be used as facet_col in plotly to split different metrics into subplots
        df["metric"] = key

    df_cum = pd.concat(dfs.values())
    # subselect rows for speed, plot has sufficient precision with 1k rows
    df_cum = df_cum.iloc[:: len(df_cum) // 1000 or 1]
    n_stable = sum(e_above_hull_true <= 0)

    if backend == "matplotlib":
        fig, axs = plt.subplots(
            ncols=min(len(metrics), 2),
            nrows=math.ceil(len(metrics) / 2),
            figsize=(15, 7),
            sharey=True,
        )
        line_kwargs = dict(
            linewidth=3, markevery=[-1], marker="x", markersize=14, markeredgewidth=2.5
        )
        # TODO breaks for len(metrics) == 1 since axs has no .flat attribute
        # assert len(axs.flat) == len(metrics), f"{len(axs.flat)=} != {len(metrics)=}"

        for metric, ax in zip(metrics, axs.flat if len(metrics) > 1 else [axs]):
            # select every n-th row of df so that 1000 rows are left for increased
            # plotting speed and reduced file size
            # falls back on every row if df has less than 1000 rows
            df = dfs[metric]
            df.iloc[:: len(df) // 1000 or 1].plot(
                ax=ax,
                legend=False,
                backend=backend,
                **line_kwargs | kwargs,
                ylabel=metric,
            )
            df = dfs[metric]
            ax.set(ylim=(0, 1), xlim=(0, None), ylabel=metric)
            for model in df_preds:
                # TODO is this if really necessary?
                if len(df[model].dropna()) == 0:
                    continue
                x_end = df[model].dropna().index[-1]
                y_end = df[model].dropna().iloc[-1]
                # add some visual guidelines to the plot
                intersect_kwargs = dict(linestyle=":", alpha=0.4, linewidth=2)
                bbox = dict(facecolor="white", alpha=0.5, edgecolor="none")
                # place model name at the end of every line
                ax.text(x_end, y_end, model, va="bottom", rotation=30, bbox=bbox)
                if "x" in project_end_point:
                    ax.plot((x_end, x_end), (0, y_end), **intersect_kwargs)
                if "y" in project_end_point:
                    ax.plot((0, x_end), (y_end, y_end), **intersect_kwargs)

            # optimal recall line finds all stable materials without any false positives
            # can be included to confirm all models achieve near optimal recall initially
            # and to see how much they overshoot n_stable
            if optimal_recall and "Recall" in metric:
                ax.plot([0, n_stable], [0, 1], color="green", linestyle="--")
                ax.text(
                    *[n_stable, 0.81],
                    optimal_recall,
                    color="green",
                    va="bottom",
                    ha="right",
                    rotation=math.degrees(math.cos(math.atan(1 / n_stable))),
                    bbox=bbox,
                )

    elif backend == "plotly":
        fig = df_cum.query(f"metric in {metrics}").plot(
            backend=backend,
            facet_col="metric",
            facet_col_wrap=2,
            facet_col_spacing=0.03,
            **kwargs,
        )

        line_kwds = dict(dash="dash", width=0.5)
        for idx, anno in enumerate(fig.layout.annotations):
            anno.text = anno.text.split("=")[1]
            anno.font.size = 16
            grid_pos = dict(row=idx // 2 + 1, col=idx % 2 + 1)
            fig.update_traces(
                hovertemplate=f"Index = %{{x:d}}<br>{anno.text} = %{{y:.2f}}",
                **grid_pos,
            )

            if optimal_recall and "recall" in anno.text.lower():
                fig.add_shape(
                    **dict(type="line", x0=0, y0=0, x1=n_stable, y1=1, **grid_pos),
                    line=line_kwds,
                )
                # annotate optimal recall line
                fig.add_annotation(
                    x=0.7 * n_stable,
                    y=0.8,
                    text=optimal_recall,
                    showarrow=False,
                    # rotate text parallel to line
                    textangle=math.degrees(math.cos(n_stable)),
                    **grid_pos,
                )
        if show_n_stable:
            fig.add_vline(x=n_stable, line=line_kwds)
            fig.add_annotation(
                x=n_stable,
                y=0.95,
                text="Stable<br>Materials",
                showarrow=False,
                xanchor="left",
                align="left",
            )
        fig.layout.legend.title = ""
        fig.update_xaxes(showticklabels=True, title="")
        fig.update_yaxes(showticklabels=True, title="")

    return fig, df_cum


def wandb_scatter(table: wandb.Table, fields: dict[str, str], **kwargs: Any) -> None:
    """Log a parity scatter plot using custom Vega spec to WandB.

    Args:
        table (wandb.Table): WandB data table.
        fields (dict[str, str]): Map from table columns to fields defined in the custom
            vega spec. Currently the only Vega fields are 'x' and 'y'.
        **kwargs: Keyword arguments passed to wandb.plot_table(string_fields=kwargs).
    """
    assert set(fields) >= {"x", "y"}, f"{fields=} must specify x and y column names"

    if "form" in fields["x"] and "form" in fields["y"]:
        kwargs.setdefault("x_label", "DFT formation energy (eV/atom)")
        kwargs.setdefault("y_label", "Predicted formation energy (eV/atom)")

    scatter_plot = wandb.plot_table(
        vega_spec_name="janosh/scatter-parity",
        data_table=table,
        fields=fields,
        string_fields=kwargs,
    )

    wandb.log({"true_pred_scatter": scatter_plot})
