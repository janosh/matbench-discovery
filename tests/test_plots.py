from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import pytest

from matbench_discovery.data import load_df_wbm_preds
from matbench_discovery.plots import (
    AxLine,
    Backend,
    WhichEnergy,
    cumulative_precision_recall,
    hist_classified_stable_vs_hull_dist,
    rolling_mae_vs_hull_dist,
)

models = ["Wrenformer", "CGCNN", "Voronoi Random Forest"]
df_wbm = load_df_wbm_preds(models, nrows=100)
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
each_pred_col = "e_above_hull_pred"
e_form_col = "e_form_per_atom_mp2020_corrected"


@pytest.mark.parametrize(
    "project_end_point,stability_threshold",
    [("", 0), ("x", 0), ("x", -0.05), ("xy", 0.1)],
)
@pytest.mark.parametrize("backend", ("matplotlib", "plotly"))
def test_cumulative_precision_recall(
    project_end_point: AxLine,
    stability_threshold: float,
    backend: Backend,
) -> None:
    fig, df_metrics = cumulative_precision_recall(
        e_above_hull_true=df_wbm[each_true_col],
        df_preds=df_wbm[models],
        backend=backend,
        project_end_point=project_end_point,
        stability_threshold=stability_threshold,
    )

    assert isinstance(df_metrics, pd.DataFrame)
    assert list(df_metrics) == models + ["metric"]

    if backend == "matplotlib":
        assert isinstance(fig, plt.Figure)
        assert all(ax.get_ylim() == (0, 1) for ax in fig.axes)
        assert (
            [ax.get_ylabel() for ax in fig.axes]
            == list(df_metrics.metric.unique())
            == ["Precision", "Recall", "F1"]
        )
    elif backend == "plotly":
        assert isinstance(fig, go.Figure)
        assert fig.layout.yaxis1.title.text == "Precision"
        assert fig.layout.yaxis2.title.text == "Recall"


@pytest.mark.parametrize("window", (0.02, 0.002))
@pytest.mark.parametrize("bin_width", (0.1, 0.001))
@pytest.mark.parametrize("x_lim", ((0, 0.6), (-0.2, 0.8)))
@pytest.mark.parametrize("backend", ("matplotlib", "plotly"))
def test_rolling_mae_vs_hull_dist(
    window: float, bin_width: float, x_lim: tuple[float, float], backend: Backend
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


@pytest.mark.parametrize("stability_threshold", (0.1, 0.01))
@pytest.mark.parametrize("x_lim", ((0, 0.6), (-0.2, 0.8)))
@pytest.mark.parametrize("which_energy", ("true", "pred"))
@pytest.mark.parametrize("backend", ("plotly", "matplotlib"))
def test_hist_classified_stable_vs_hull_dist(
    stability_threshold: float,
    x_lim: tuple[float, float],
    which_energy: WhichEnergy,
    backend: Backend,
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    df_wbm[each_pred_col] = (
        df_wbm[each_true_col] + df_wbm[models[0]] - df_wbm[e_form_col]
    )

    ax = hist_classified_stable_vs_hull_dist(
        df_wbm,
        each_true_col=each_true_col,
        each_pred_col=each_pred_col,
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
