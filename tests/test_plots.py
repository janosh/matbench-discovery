from typing import Literal
from unittest.mock import MagicMock, patch

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
    wandb_scatter,
)

AxLine = Literal["x", "y", "xy", ""]
models = ["MEGNet", "CGCNN", "Voronoi RF"]
df_wbm = load_df_wbm_with_preds(models=models, nrows=100)
# TODO remove pd.DataFrame type cast pending https://github.com/astral-sh/ty/issues/1075
df_preds = pd.DataFrame(df_wbm[models])


@pytest.mark.parametrize(
    "stability_threshold, metrics",
    [
        (0, ("Recall",)),
        (-0.05, ("Recall", "MAE")),
        (0.1, ("Recall", "Precision", "RMSE")),
    ],
)
def test_cumulative_metrics(
    stability_threshold: float,
    metrics: tuple[str, ...],
) -> None:
    fig, df_cumu_metrics = cumulative_metrics(
        e_above_hull_true=df_wbm[MbdKey.each_true],
        df_preds=df_preds,
        stability_threshold=stability_threshold,
        metrics=metrics,
    )

    assert isinstance(df_cumu_metrics, pd.DataFrame)
    assert list(df_cumu_metrics) == ["index", *models, "metric"]

    assert isinstance(fig, go.Figure)
    subplot_titles = {anno.text.split("=")[-1] for anno in fig.layout.annotations}
    assert subplot_titles >= set(metrics)


def test_cumulative_metrics_raises() -> None:
    with pytest.raises(
        ValueError,
        match=r"invalid_metrics=\{'invalid'\}, should be case-insensitive subset of",
    ):
        cumulative_metrics(
            e_above_hull_true=df_wbm[MbdKey.each_true],
            df_preds=df_preds,
            metrics=("invalid",),
        )


@pytest.mark.parametrize(
    "window, bin_width, x_lim, show_dft_acc",
    [
        (0.02, 0.1, (0, 0.6), True),
        (0.002, 0.001, (-0.2, 0.8), False),
    ],
)
def test_rolling_mae_vs_hull_dist(
    window: float,
    bin_width: float,
    x_lim: tuple[float, float],
    show_dft_acc: bool,
) -> None:
    ax, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_wbm[models[0]],
        e_above_hull_preds=df_preds,
        x_lim=x_lim,
        show_dft_acc=show_dft_acc,
        window=window,
        bin_width=bin_width,
    )

    assert isinstance(df_err, pd.DataFrame)
    assert isinstance(df_std, pd.DataFrame)
    assert list(df_err) == list(df_std) == models, (
        f"expected {list(df_err)} == {list(df_std)} == {models}"
    )

    expected_ylabel = "Rolling MAE (eV/atom)"
    assert isinstance(ax, go.Figure)
    assert ax.layout.yaxis.title.text == expected_ylabel
    assert ax.layout.xaxis.title.text == "E<sub>above MP hull</sub> (eV/atom)"


@pytest.mark.parametrize(
    "stability_threshold,x_lim,which_energy",
    [
        (0.1, (0, 0.6), "true"),
        (0.01, (-0.2, 0.8), "pred"),
        (0.1, (-0.2, 0.8), "true"),
        (0.01, (0, 0.6), "pred"),
    ],
)
def test_hist_classified_stable_vs_hull_dist(
    stability_threshold: float,
    x_lim: tuple[float, float],
    which_energy: Literal["true", "pred"],
) -> None:
    df_wbm[Key.each_pred] = (
        df_wbm[MbdKey.each_true] + df_wbm[models[0]] - df_wbm[MbdKey.e_form_dft]
    )

    fig = hist_classified_stable_vs_hull_dist(
        df_wbm,
        each_true_col=MbdKey.each_true,
        each_pred_col=Key.each_pred,
        stability_threshold=stability_threshold,
        x_lim=x_lim,
        which_energy=which_energy,
    )

    assert isinstance(fig, go.Figure)
    assert fig.layout.yaxis.title.text == "Count"


def test_plotly_markers_line_styles() -> None:
    """Test plotly markers and line styles are populated correctly."""
    assert len(plotly_markers) > 100
    assert len(plotly_line_styles) > 100
    assert {*map(type, plotly_markers)} == {str}, "expect all markers are strings"
    assert {*map(type, plotly_line_styles)} == {str}
    assert "longdashdot" in plotly_line_styles
    assert "circle" in plotly_markers


@pytest.mark.parametrize(
    "legend_loc, should_raise",
    [("figure", False), ("below", False), ("default", False), ("invalid", True)],
)
def test_rolling_mae_vs_hull_dist_legend_loc(
    legend_loc: str, should_raise: bool
) -> None:
    """Test rolling_mae_vs_hull_dist with valid/invalid legend locations."""
    call_fn = lambda: rolling_mae_vs_hull_dist(
        e_above_hull_true=df_wbm[models[0]],
        e_above_hull_preds=df_preds,
        x_lim=(0, 0.6),
        window=0.02,
        bin_width=0.1,
        legend_loc=legend_loc,  # type: ignore[arg-type]
    )
    if should_raise:
        with pytest.raises(ValueError, match=f"Unexpected legend_loc='{legend_loc}'"):
            call_fn()
    else:
        fig, df_err, df_std = call_fn()
        assert isinstance(fig, go.Figure)
        assert isinstance(df_err, pd.DataFrame)
        assert isinstance(df_std, pd.DataFrame)


@pytest.mark.parametrize(
    "fields, should_raise, check_labels",
    [
        ({"x": "col_x"}, True, False),  # missing 'y'
        ({"y": "col_y"}, True, False),  # missing 'x'
        ({"x": "col_x", "y": "col_y"}, False, False),  # basic valid case
        ({"x": "e_form_true", "y": "e_form_pred"}, False, True),  # formation energy
    ],
)
def test_wandb_scatter(fields: dict, should_raise: bool, check_labels: bool) -> None:
    """Test wandb_scatter with valid/invalid fields and label defaults."""
    mock_table = MagicMock()

    if should_raise:
        with pytest.raises(ValueError, match="must specify x=str and y=str"):
            wandb_scatter(mock_table, fields)
    else:
        with patch("matbench_discovery.plots.wandb") as mock_wandb:
            mock_wandb.plot_table.return_value = MagicMock()
            wandb_scatter(mock_table, fields)

            mock_wandb.plot_table.assert_called_once()
            mock_wandb.log.assert_called_once()

            if check_labels:
                call_kwargs = mock_wandb.plot_table.call_args.kwargs
                assert call_kwargs["string_fields"]["x_label"] == (
                    "DFT formation energy (eV/atom)"
                )
                assert call_kwargs["string_fields"]["y_label"] == (
                    "Predicted formation energy (eV/atom)"
                )
