from typing import Literal

import pandas as pd
import plotly.graph_objects as go
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey
from matbench_discovery.plots import (
    cumulative_metrics,
    hist_classified_stable_vs_hull_dist,
    plotly_line_styles,
    plotly_markers,
    rolling_mae_vs_hull_dist,
)

AxLine = Literal["x", "y", "xy", ""]
models = ["MEGNet", "CGCNN", "Voronoi RF"]
df_wbm = load_df_wbm_with_preds(models=models, nrows=100)


@pytest.mark.parametrize(
    "stability_threshold",
    [(0), (-0.05), (0.1)],
)
@pytest.mark.parametrize(
    "metrics",
    [("Recall",), ("Recall", "MAE"), ("Recall", "Precision", "RMSE")],
)
def test_cumulative_metrics(
    stability_threshold: float,
    metrics: tuple[str, ...],
) -> None:
    fig, df_cumu_metrics = cumulative_metrics(
        e_above_hull_true=df_wbm[MbdKey.each_true],
        df_preds=df_wbm[models],
        stability_threshold=stability_threshold,
        metrics=metrics,
    )

    assert isinstance(df_cumu_metrics, pd.DataFrame)
    assert list(df_cumu_metrics) == [*models, "metric"]

    assert isinstance(fig, go.Figure)
    subplot_titles = {anno.text.split("=")[-1] for anno in fig.layout.annotations}
    assert subplot_titles >= set(metrics)


def test_cumulative_metrics_raises() -> None:
    with pytest.raises(
        ValueError,
        match="invalid_metrics={'invalid'}, should be case-insensitive subset of",
    ):
        cumulative_metrics(
            e_above_hull_true=df_wbm[MbdKey.each_true],
            df_preds=df_wbm[models],
            metrics=("invalid",),
        )


@pytest.mark.parametrize("window", [0.02, 0.002])
@pytest.mark.parametrize("bin_width", [0.1, 0.001])
@pytest.mark.parametrize("x_lim", [(0, 0.6), (-0.2, 0.8)])
@pytest.mark.parametrize("show_dft_acc", [True, False])
def test_rolling_mae_vs_hull_dist(
    window: float,
    bin_width: float,
    x_lim: tuple[float, float],
    show_dft_acc: bool,
) -> None:
    ax, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_wbm[models[0]],
        e_above_hull_preds=df_wbm[models],
        x_lim=x_lim,
        show_dft_acc=show_dft_acc,
        window=window,
        bin_width=bin_width,
    )

    assert isinstance(df_err, pd.DataFrame)
    assert isinstance(df_std, pd.DataFrame)
    assert (
        list(df_err) == list(df_std) == models
    ), f"expected {list(df_err)} == {list(df_std)} == {models}"

    expected_ylabel = "Rolling MAE (eV/atom)"
    assert isinstance(ax, go.Figure)
    assert ax.layout.yaxis.title.text == expected_ylabel
    assert ax.layout.xaxis.title.text == "E<sub>above MP hull</sub> (eV/atom)"


@pytest.mark.parametrize("stability_threshold", [0.1, 0.01])
@pytest.mark.parametrize("x_lim", [(0, 0.6), (-0.2, 0.8)])
@pytest.mark.parametrize("which_energy", ["true", "pred"])
def test_hist_classified_stable_vs_hull_dist(
    stability_threshold: float,
    x_lim: tuple[float, float],
    which_energy: Literal["true", "pred"],
) -> None:
    df_wbm[Key.each_pred] = (
        df_wbm[MbdKey.each_true] + df_wbm[models[0]] - df_wbm[MbdKey.e_form_dft]
    )

    ax = hist_classified_stable_vs_hull_dist(
        df_wbm,
        each_true_col=MbdKey.each_true,
        each_pred_col=Key.each_pred,
        stability_threshold=stability_threshold,
        x_lim=x_lim,
        which_energy=which_energy,
    )

    assert isinstance(ax, go.Figure)
    assert ax.layout.yaxis.title.text == "count"


def test_plotly_markers_line_styles() -> None:
    assert len(plotly_markers) > 100
    assert len(plotly_line_styles) > 100
    assert {*map(type, plotly_markers)} == {str}, "expect all markers are strings"
    assert {*map(type, plotly_line_styles)} == {str}
    assert "longdashdot" in plotly_line_styles
    assert "circle" in plotly_markers


if __name__ == "__main__":
    pytest.main(["-vv", __file__])
