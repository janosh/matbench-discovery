"""Plotting functions for analyzing model performance on materials discovery."""

from __future__ import annotations

import functools
import math
from collections import defaultdict
from typing import TYPE_CHECKING, Any, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import scipy.interpolate
import scipy.stats
import wandb
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from plotly.validators.scatter.line import DashValidator
from plotly.validators.scatter.marker import SymbolValidator
from tqdm import tqdm

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.metrics import classify_stable

if TYPE_CHECKING:
    from collections.abc import Sequence

__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"

Backend = Literal["matplotlib", "plotly"]

plotly_markers = SymbolValidator().values[2::3]  # noqa: PD011
plotly_line_styles = DashValidator().values[:-1]  # noqa: PD011
plotly_colors = px.colors.qualitative.Plotly
# repeat line styles/colors as many as times as needed to match number of markers
plotly_line_styles *= len(plotly_markers) // len(plotly_line_styles)
plotly_colors *= len(plotly_markers) // len(plotly_colors)


# color list https://plotly.com/python-api-reference/generated/plotly.graph_objects.layout
colorway = ("lightseagreen", "orange", "lightsalmon", "dodgerblue")
clf_labels = ("True Positive", "False Negative", "False Positive", "True Negative")
clf_colors = ("lightseagreen", "orange", "lightsalmon", "dodgerblue")
clf_color_map = dict(zip(clf_labels, clf_colors))


def hist_classified_stable_vs_hull_dist(
    df: pd.DataFrame,
    each_true_col: str,
    each_pred_col: str,
    ax: plt.Axes | None = None,
    which_energy: Literal["true", "pred"] = "true",
    stability_threshold: float | None = 0,
    x_lim: tuple[float, float] = (-0.7, 0.7),
    n_bins: int = 200,
    rolling_acc: float | None = 0.02,
    backend: Backend = "plotly",
    y_label: str = "Number of materials",
    clf_labels: Sequence[str] = clf_labels,
    **kwargs: Any,
) -> plt.Axes | go.Figure:
    """Histogram of the energy difference (either according to DFT ground truth - the
    default - or the model predicted energy) to the convex hull for materials in the
    WBM data set. The histogram is broken down into true positives, false negatives,
    false positives, and true negatives based on whether the model predicts candidates
    to be below the known convex hull. Ideally, in discovery setting a model should
    exhibit high recall, i.e. the majority of materials below the convex hull being
    correctly identified by the model.

    See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.

    Args:
        df (pd.DataFrame): Data frame containing true and predicted hull distances.
        each_true_col (str): Name of column with energy above convex hull according to
            DFT ground truth (in eV / atom).
        each_pred_col (str): Name of column with energy above convex hull predicted by
            model (in eV / atom). Same as true energy to convex hull plus predicted
            minus true formation energy.
        ax (plt.Axes, optional): matplotlib axes to plot on.
        which_energy ('true' | 'pred', optional): Whether to use the true (DFT) hull
            distance or the model's predicted hull distance for the histogram.
        stability_threshold (float, optional): set stability threshold as distance to
            convex hull in eV/atom, usually 0 or 0.1 eV.
        x_lim (tuple[float, float]): x-axis limits.
        n_bins (int): Number of bins in histogram.
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
    x_col = dict(true=each_true_col, pred=each_pred_col)[which_energy]

    df_plot = pd.DataFrame()
    for facet, df_group in (
        df.groupby(kwargs["facet_col"]) if "facet_col" in kwargs else [(None, df)]
    ):
        true_pos, false_neg, false_pos, true_neg = classify_stable(
            df_group[each_true_col], df_group[each_pred_col], stability_threshold
        )

        # switch between hist of DFT-computed and model-predicted convex hull distance
        e_above_hull = df_group[x_col]
        each_true_pos = e_above_hull[true_pos]
        each_true_neg = e_above_hull[true_neg]
        each_false_neg = e_above_hull[false_neg]
        each_false_pos = e_above_hull[false_pos]
        # n_true_pos, n_false_pos, n_true_neg, n_false_neg = map(
        #     sum, (true_pos, false_pos, true_neg, false_neg)
        # )

        df_group[(clf_col := "classified")] = np.array(clf_labels)[
            true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
        ]

        # calculate histograms for each category
        hist_true_pos, bin_edges = np.histogram(
            each_true_pos - 0.001, bins=n_bins, range=x_lim
        )
        hist_true_neg, _ = np.histogram(each_true_neg + 0.001, bins=n_bins, range=x_lim)
        hist_false_neg, _ = np.histogram(
            each_false_neg - 0.001, bins=n_bins, range=x_lim
        )
        hist_false_pos, _ = np.histogram(
            each_false_pos + 0.001, bins=n_bins, range=x_lim
        )

        # combine histograms into a single dataframe
        df_hist = pd.DataFrame(
            (hist_true_pos, hist_false_neg, hist_false_pos, hist_true_neg),
            index=clf_labels,
        ).T
        df_hist[x_col] = bin_edges[:-1]
        value_name = "count"
        df_melt = df_hist.melt(
            id_vars=x_col,
            value_vars=clf_labels,
            var_name=clf_col,
            value_name=value_name,
        ).assign(**{kwargs.get("facet_col", ""): facet})
        df_plot = pd.concat([df_plot, df_melt])

    if backend == "plotly":
        kwargs.update(
            barmode="stack",
            range_x=x_lim,
            color_discrete_map=clf_color_map,
            x=x_col,
            color=clf_col,
            y=value_name,
        )
    else:
        kwargs.update(
            column=[x_col],
            legend=False,
            ax=ax,
            xlim=x_lim,
            # color=df[clf_col],
        )

    fig = df_plot.round(4).plot.bar(backend=backend, **kwargs)

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
        fig.layout.legend.update(title=None, y=0.5, xanchor="right", x=1)
        fig.update_traces(marker_line=dict(color="rgba(0, 0, 0, 0)"))
        fig.update_layout(bargap=0)
        fig.update_xaxes(matches=None)
        fig.update_yaxes(matches=None)

        if stability_threshold is not None:
            anno = dict(text="Stability threshold", xanchor="center", yanchor="bottom")
            fig.add_vline(stability_threshold, line=dict(dash="dash"), annotation=anno)

    if rolling_acc:
        # add moving average of the accuracy computed within given window
        # as a function of e_above_hull shown as blue line (right axis)

        # --- moving average of the accuracy
        # compute rolling accuracy in rolling_acc-sized intervals
        bins = np.arange(df[x_col].min(), df[x_col].max(), rolling_acc)
        bin_counts = np.histogram(df[each_true_col], bins)[0]
        bin_true_pos = np.histogram(each_true_pos, bins)[0]
        bin_true_neg = np.histogram(each_true_neg, bins)[0]

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
    df_rolling_err: pd.DataFrame | None = None,
    df_err_std: pd.DataFrame | None = None,
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
    pbar: bool = True,
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
        show_dft_acc (bool, optional): If True, change color of the triangle of peril's
            tip and annotate it with 'Corrected GGA Accuracy' at rolling MAE of 25
            meV/atom. Defaults to False.
        show_dummy_mae (bool, optional): If True, plot a line at the dummy MAE of always
            predicting the target mean.
        pbar (bool, optional): If True, show a progress bar during rolling MAE
            calculation. Defaults to True.
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

        for model in (
            prog_bar := tqdm(models, desc="Calculating rolling MAE", disable=not pbar)
        ):
            prog_bar.set_postfix_str(model)
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

    fig = df_rolling_err.plot(backend=backend, **kwargs)

    if just_plot_lines:
        # return earlier if all plot objects besides the line were already drawn by a
        # previous call
        return fig, df_rolling_err, df_err_std

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
            fig.transData,
            window,
            window_bar_anno,
            "lower left",
            pad=0.5,
            frameon=False,
            size_vertical=0.002,
        )
        # indicate size of MAE averaging window
        fig.add_artist(scale_bar)

        fig.fill_between(
            (-1, -dft_acc, dft_acc, 1) if show_dft_acc else (-1, 0, 1),
            (1, 1, 1, 1) if show_dft_acc else (1, 1, 1),
            (1, dft_acc, dft_acc, 1) if show_dft_acc else (1, 0, 1),
            color="tab:red",
            alpha=0.2,
        )

        if show_dft_acc:
            fig.fill_between(
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
            fig.annotate(
                xy=(-dft_acc, dft_acc),
                xytext=(-2 * dft_acc, dft_acc),
                text="Corrected GGA\nAccuracy",
                arrowprops=arrowprops,
                verticalalignment="center",
                horizontalalignment="right",
            )

        fig.axhline(dummy_mae, color="tab:blue", linestyle="--", linewidth=0.5)
        fig.text(dummy_mae, 0.1, dummy_mae_text)

        fig.text(
            0, 0.13, r"MAE > $|E_\mathrm{above\ hull}|$", horizontalalignment="center"
        )
        fig.set(xlabel=r"$E_\mathrm{above\ hull}$ (eV/atom)", ylabel=y_label)
        fig.set(xlim=x_lim, ylim=y_lim)
        plt_line_styles = "- -- -. :".split() * 10
        plt_markers = "o s ^ v D * p X".split() * 10

        for idx, line in enumerate(fig.lines):
            line_label = line.get_label()
            if line_label.startswith("_"):
                continue
            ls, marker = plt_line_styles[idx], plt_markers[idx]
            line.set(ls=ls, marker=marker, markeredgewidth=0.5, markeredgecolor="black")
            line.set_markevery(8)

    elif backend == "plotly":
        for idx, model in enumerate(df_rolling_err if with_sem else []):
            # set legendgroup to model name so SEM shading toggles with model curve
            fig.data[idx].legendgroup = model
            # set SEM area to same color as model curve
            fig.add_scatter(
                x=list(bins) + list(bins)[::-1],  # bins, then bins reversed
                # upper, then lower reversed
                y=list(df_rolling_err[model] + 3 * df_err_std[model])
                + list(df_rolling_err[model] - 3 * df_err_std[model])[::-1],
                mode="lines",
                line=dict(color="white", width=0),
                fill="toself",
                legendgroup=model,
                fillcolor=fig.data[idx].line.color,
                opacity=0.4,
                showlegend=False,
            )

        fig.layout.legend.update(title="", x=0, y=0, yanchor="bottom")
        # change tooltip precision to 2 decimal places
        fig.update_traces(hovertemplate="x = %{x:.2f} eV/atom<br>y = %{y:.2f} eV/atom")
        fig.layout.xaxis.title.text = "E<sub>above MP hull</sub> (eV/atom)"
        fig.layout.yaxis.title.text = "rolling MAE (eV/atom)"
        fig.update_xaxes(range=x_lim)
        fig.update_yaxes(range=y_lim)
        # exclude from hover tooltip
        scatter_kwds = dict(
            fill="toself", opacity=0.2, hoverinfo="skip", showlegend=False
        )
        triangle_anno = "MAE > |E<sub>hull dist</sub>|"
        fig.add_scatter(
            x=(-1, -dft_acc, dft_acc, 1) if show_dft_acc else (-1, 0, 1),
            y=(1, dft_acc, dft_acc, 1) if show_dft_acc else (1, 0, 1),
            name=triangle_anno,
            fillcolor="red",
            # remove triangle border
            line=dict(color="rgba(0,0,0,0)"),
            **scatter_kwds,
        )
        fig.add_annotation(
            x=0, y=0.7, text=triangle_anno, showarrow=False, yref="paper"
        )

        if show_dummy_mae:
            fig.add_hline(
                y=dummy_mae,
                line=dict(dash="dash", width=0.5),
                annotation_text=dummy_mae_text,
            )

        if show_dft_acc:
            fig.add_scatter(
                x=(-dft_acc, dft_acc, 0, -dft_acc),
                y=(dft_acc, dft_acc, 0, dft_acc),
                name="MAE < |Corrected GGA error|",
                fillcolor="red",
                **scatter_kwds,
            )
            fig.add_annotation(
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

        fig.data = fig.data[::-1]  # bring px.line() to front
        # plot rectangle to indicate MAE window size
        x0, y0 = x_lim[1] - 0.01, y_lim[0] + 0.01
        fig.add_annotation(
            x=x0 - window,
            y=y0,
            text=window_bar_anno,
            showarrow=False,
            xshift=-4,
            yshift=-6,
            yanchor="bottom",
            xanchor="right",
        )
        fig.add_shape(type="rect", x0=x0, y0=y0, x1=x0 - window, y1=y0 + window / 5)

        from matbench_discovery.preds import model_styles

        for trace in fig.data:
            if style := model_styles.get(trace.name):
                ls, _marker, color = style
                trace.line = dict(color=color, dash=ls, width=2)
            # marker_spacing = 2
            # fig.add_scatter(
            #     x=trace.x[::marker_spacing],
            #     y=trace.y[::marker_spacing],
            #     mode="markers",
            #     marker=dict(symbol=marker, color=trace.line.color, size=8),
            #     showlegend=False,
            #     legendgroup=getattr(trace, "legendgroup", None),
            # )
    return fig, df_rolling_err, df_err_std


def cumulative_metrics(
    e_above_hull_true: pd.Series,
    df_preds: pd.DataFrame,
    metrics: Sequence[str] = ("Precision", "Recall"),
    stability_threshold: float = 0,  # set stability threshold as distance to convex
    # hull in eV / atom, usually 0 or 0.1 eV
    project_end_point: Literal["x", "y", "xy", ""] = "xy",
    optimal_recall: str | None = "Optimal Recall",
    show_n_stable: bool = True,
    backend: Backend = "plotly",
    n_points: int = 100,
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
        metrics (Sequence[str], optional): Which metrics to plot. Any subset of
            ("Precision", "Recall", "F1", "MAE", "RMSE").
            Defaults to ('Precision', 'Recall').
        stability_threshold (float, optional): Max distance above convex hull before
            material is considered unstable. Defaults to 0.
        project_end_point ('x' | 'y' | 'xy' | '', optional): Whether to project end
        points of precision and recall curves to the x/y axis. Defaults to '', i.e. no
            axis projection lines.
        optimal_recall (str | None, optional): Label for the optimal recall line.
            Defaults to 'Optimal Recall'. Set to None to not plot the line.
        show_n_stable (bool, optional): Whether to show a horizontal line at the true
            number of stable materials. Defaults to True.
        backend ('matplotlib' | 'plotly'], optional): Which plotting engine to use.
            Changes the return type. Defaults to 'plotly'.
        n_points (int, optional): Number of points to use for interpolation of the
            metric curves. Defaults to 100.
        **kwargs: Keyword arguments passed to df.plot().

    Returns:
        tuple[plt.Figure | go.Figure, pd.DataFrame]: The matplotlib/plotly figure and
            dataframe of cumulative metrics for each model.
    """
    dfs: dict[str, pd.DataFrame] = defaultdict(pd.DataFrame)

    # largest number of materials predicted stable by any model, determines x-axis range
    n_max_pred_stable = (df_preds < stability_threshold).sum().max()
    # use log2-spaced sampling to get higher sampling density at equal file size for
    # start of the discovery campaign where model performance fluctuates more
    longest_xs = np.logspace(0, np.log2(n_max_pred_stable - 1), n_points, base=2)
    for metric in metrics:
        dfs[metric].index = longest_xs

    valid_metrics = {"Precision", "Recall", "F1", "MAE", "RMSE"}
    if invalid_metrics := set(metrics) - valid_metrics:
        raise ValueError(
            f"{invalid_metrics=}, should be case-insensitive subset of {valid_metrics=}"
        )

    for model_name in df_preds:
        each_pred = df_preds[model_name].sort_values()
        # sort targets by model ranking
        each_true = e_above_hull_true.loc[each_pred.index]

        true_pos_cum, false_neg_cum, false_pos_cum, true_neg_cum = map(
            np.cumsum, classify_stable(each_true, each_pred, stability_threshold)
        )

        # precision aka positive predictive value (PPV)
        n_total_pos = true_pos_cum.iloc[-1] + false_neg_cum.iloc[-1]
        precision_cum = true_pos_cum / (true_pos_cum + false_pos_cum)
        recall_cum = true_pos_cum / n_total_pos  # aka true_pos_rate aka sensitivity

        n_pred_stable = sum(each_pred <= stability_threshold)
        model_range = np.arange(n_pred_stable)  # xs for interpolation
        xs_model = longest_xs[longest_xs < n_pred_stable - 1]  # xs for plotting

        cubic_interpolate = functools.partial(scipy.interpolate.interp1d, kind="cubic")

        if "Precision" in metrics:
            prec_interp = cubic_interpolate(model_range, precision_cum[:n_pred_stable])
            dfs["Precision"][model_name] = dict(zip(xs_model, prec_interp(xs_model)))

        if "Recall" in metrics:
            recall_interp = cubic_interpolate(model_range, recall_cum[:n_pred_stable])
            dfs["Recall"][model_name] = dict(zip(xs_model, recall_interp(xs_model)))

        if "F1" in metrics:
            f1_cum = 2 * (precision_cum * recall_cum) / (precision_cum + recall_cum)
            f1_interp = cubic_interpolate(model_range, f1_cum[:n_pred_stable])
            dfs["F1"][model_name] = dict(zip(xs_model, f1_interp(xs_model)))

        cum_counts = np.arange(1, len(each_true) + 1)
        if "MAE" in metrics:
            cum_errors = (each_true - each_pred).abs().cumsum()
            mae_cum = cum_errors / cum_counts
            mae_interp = cubic_interpolate(model_range, mae_cum[:n_pred_stable])
            dfs["MAE"][model_name] = dict(zip(xs_model, mae_interp(xs_model)))

        if "RMSE" in metrics:
            rmse_cum = (((each_true - each_pred) ** 2).cumsum() / cum_counts) ** 0.5
            rmse_interp = cubic_interpolate(model_range, rmse_cum[:n_pred_stable])
            dfs["RMSE"][model_name] = dict(zip(xs_model, rmse_interp(xs_model)))

    for key in dfs:
        # drop all-NaN rows so plotly plot x-axis only extends to largest number of
        # predicted materials by any model
        dfs[key] = dfs[key].dropna(how="all")
        # will be used as facet_col in plotly to split different metrics into subplots
        dfs[key]["metric"] = key

    df_cum = pd.concat(dfs.values())
    # subselect rows for speed, plot has sufficient precision with 1k rows
    n_stable = sum(e_above_hull_true <= STABILITY_THRESHOLD)

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
            df.plot(
                ax=ax,
                legend=False,
                backend=backend,
                **line_kwargs | kwargs,
                ylabel=metric,
            )
            df = dfs[metric]
            ax.set(ylim=(0, 1), xlim=(0, None), ylabel=metric)
            for model in df_preds:
                # TODO is this really necessary?
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
            # can be included to confirm all models achieve near optimal recall
            # initially and to see how much they overshoot n_stable
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
        n_cols = kwargs.pop("facet_col_wrap", 2)
        kwargs.setdefault("facet_col_spacing", 0.03)
        fig = df_cum.plot(
            backend=backend,
            facet_col="metric",
            facet_col_wrap=n_cols,
            **kwargs,
        )

        line_kwds = dict(dash="dash", width=0.5)
        for idx, anno in enumerate(fig.layout.annotations):
            anno.text = anno.text.split("=")[1]
            anno.font.size = 16
            grid_pos = dict(row=idx // n_cols + 1, col=idx % n_cols + 1)
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
                    # angle not quite right, could be improved
                    textangle=math.degrees(math.cos(n_stable)),
                    **grid_pos,
                )
        if show_n_stable:
            fig.add_vline(x=n_stable, line=line_kwds)
            fig.add_annotation(
                x=n_stable,
                y=0.95,
                text="Stable Materials",
                showarrow=False,
                xanchor="left",
                align="left",
            )
        fig.layout.legend.title = ""
        fig.update_xaxes(showticklabels=True, title="", matches=None)
        fig.update_yaxes(showticklabels=True, title="", matches=None)

    return fig, df_cum


def wandb_scatter(table: wandb.Table, fields: dict[str, str], **kwargs: Any) -> None:
    """Log a parity scatter plot using custom Vega spec to WandB.

    Args:
        table (wandb.Table): WandB data table.
        fields (dict[str, str]): Map from table columns to fields defined in the custom
            vega spec. Currently the only Vega fields are 'x' and 'y'.
        **kwargs: Keyword arguments passed to wandb.plot_table(string_fields=kwargs).
    """
    if set(fields) < {"x", "y"}:
        raise ValueError(f"{fields=} must specify x=str and y=str column names")

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
