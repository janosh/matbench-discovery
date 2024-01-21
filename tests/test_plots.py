from __future__ import annotations

from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import pytest

from matbench_discovery import Key
from matbench_discovery.plots import (
    Backend,
    cumulative_metrics,
    hist_classified_stable_vs_hull_dist,
    plotly_line_styles,
    plotly_markers,
    rolling_mae_vs_hull_dist,
)
from matbench_discovery.preds import load_df_wbm_with_preds

AxLine = Literal["x", "y", "xy", ""]
models = ["MEGNet", "CGCNN", "Voronoi RF"]
df_wbm = load_df_wbm_with_preds(models, nrows=100)


@pytest.mark.parametrize(
    "project_end_point,stability_threshold",
    [("", 0), ("x", 0), ("x", -0.05), ("xy", 0.1)],
)
@pytest.mark.parametrize("backend", ["matplotlib", "plotly"])
@pytest.mark.parametrize(
    "metrics",
    [("Recall",), ("Recall", "MAE"), ("Recall", "Precision", "RMSE")],
)
def test_cumulative_metrics(
    project_end_point: AxLine,
    stability_threshold: float,
    backend: Backend,
    metrics: tuple[str, ...],
) -> None:
    fig, df_metrics = cumulative_metrics(
        e_above_hull_true=df_wbm[Key.each_true],
        df_preds=df_wbm[models],
        backend=backend,
        project_end_point=project_end_point,
        stability_threshold=stability_threshold,
        metrics=metrics,
    )

    assert isinstance(df_metrics, pd.DataFrame)
    assert list(df_metrics) == [*models, "metric"]

    if backend == "matplotlib":
        assert isinstance(fig, plt.Figure)
        assert all(ax.get_ylim() == (0, 1) for ax in fig.axes)
        assert {ax.get_ylabel() for ax in fig.axes} >= {*metrics}
    elif backend == "plotly":
        assert isinstance(fig, go.Figure)
        subplot_titles = {anno.text.split("=")[-1] for anno in fig.layout.annotations}
        assert subplot_titles >= set(metrics)


def test_cumulative_metrics_raises() -> None:
    with pytest.raises(
        ValueError,
        match="invalid_metrics={'invalid'}, should be case-insensitive subset of",
    ):
        cumulative_metrics(
            e_above_hull_true=df_wbm[Key.each_true],
            df_preds=df_wbm[models],
            metrics=("invalid",),
        )


@pytest.mark.parametrize("window", [0.02, 0.002])
@pytest.mark.parametrize("bin_width", [0.1, 0.001])
@pytest.mark.parametrize("x_lim", [(0, 0.6), (-0.2, 0.8)])
@pytest.mark.parametrize("backend", ["matplotlib", "plotly"])
@pytest.mark.parametrize("show_dft_acc", [True, False])
def test_rolling_mae_vs_hull_dist(
    window: float,
    bin_width: float,
    x_lim: tuple[float, float],
    backend: Backend,
    show_dft_acc: bool,
) -> None:
    kwargs = dict(window=window, bin_width=bin_width, backend=backend)
    if backend == "matplotlib":
        ax = plt.figure().gca()  # new figure ensures test functions use different axes
        kwargs["ax"] = ax

    for model_name in models:
        ax, df_err, df_std = rolling_mae_vs_hull_dist(
            e_above_hull_true=df_wbm[model_name],
            e_above_hull_errors=df_wbm[models],
            x_lim=x_lim,
            show_dft_acc=show_dft_acc,
            **kwargs,  # type: ignore[arg-type]
        )

    assert isinstance(df_err, pd.DataFrame)
    assert isinstance(df_std, pd.DataFrame)
    assert (
        list(df_err) == list(df_std) == models
    ), f"expected {list(df_err)} == {list(df_std)} == {models}"

    expected_ylabel = "rolling MAE (eV/atom)"
    if backend == "matplotlib":
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylim()[0] == 0
        assert ax.get_ylabel() == expected_ylabel
        assert ax.get_xlabel() == r"$E_\mathrm{above\ hull}$ (eV/atom)"
    elif backend == "plotly":
        assert isinstance(ax, go.Figure)
        assert ax.layout.yaxis.title.text == expected_ylabel
        assert ax.layout.xaxis.title.text == "E<sub>above MP hull</sub> (eV/atom)"


@pytest.mark.parametrize("stability_threshold", [0.1, 0.01])
@pytest.mark.parametrize("x_lim", [(0, 0.6), (-0.2, 0.8)])
@pytest.mark.parametrize("which_energy", ["true", "pred"])
@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
def test_hist_classified_stable_vs_hull_dist(
    stability_threshold: float,
    x_lim: tuple[float, float],
    which_energy: Literal["true", "pred"],
    backend: Backend,
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    df_wbm[Key.each_pred] = (
        df_wbm[Key.each_true] + df_wbm[models[0]] - df_wbm[Key.e_form]
    )

    ax = hist_classified_stable_vs_hull_dist(
        df_wbm,
        each_true_col=Key.each_true,
        each_pred_col=Key.each_pred,
        ax=ax,
        stability_threshold=stability_threshold,
        x_lim=x_lim,
        which_energy=which_energy,
        backend=backend,
    )

    if backend == "matplotlib":
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylabel() == "Number of materials"
    else:
        assert isinstance(ax, go.Figure)
        assert ax.layout.yaxis.title.text == "count"


def test_plotly_markers_line_styles() -> None:
    assert len(plotly_markers) > 100
    assert len(plotly_line_styles) > 100
    assert {*map(type, plotly_markers)} == {str}, "expect all markers are strings"
    assert {*map(type, plotly_line_styles)} == {str}
    assert "longdashdot" in plotly_line_styles
    assert "circle" in plotly_markers
