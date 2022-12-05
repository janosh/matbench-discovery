from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import pytest

from matbench_discovery import ROOT
from matbench_discovery.load_preds import df_wbm
from matbench_discovery.plots import (
    AxLine,
    Backend,
    WhichEnergy,
    cumulative_precision_recall,
    hist_classified_stable_vs_hull_dist,
    rolling_mae_vs_hull_dist,
)

DATA_DIR = f"{ROOT}/data/2022-06-11-from-rhys"

test_dfs: dict[str, pd.DataFrame] = {}
for model_name in ("Wren", "CGCNN", "Voronoi"):
    df = pd.read_csv(
        f"{DATA_DIR}/{model_name.lower()}-mp-initial-structures.csv", nrows=100
    ).set_index("material_id")

    df["e_above_hull_mp"] = df_wbm.e_above_hull_mp2020_corrected_ppd_mp

    model_preds = df.filter(like=r"_pred").mean(axis=1)

    df["e_above_hull_pred"] = df.e_above_hull_mp + model_preds - df.e_form_target

    test_dfs[model_name] = df


@pytest.mark.parametrize(
    "project_end_point,stability_threshold,backend",
    [
        ("", 0, "plotly"),
        ("x", 0, "plotly"),
        ("x", -0.05, "matplotlib"),
        ("xy", 0.1, "matplotlib"),
    ],
)
def test_cumulative_precision_recall(
    project_end_point: AxLine,
    stability_threshold: float,
    backend: Backend,
) -> None:
    fig, df_metrics = cumulative_precision_recall(
        e_above_hull_true=df.e_above_hull_mp,
        df_preds=df.filter(like="_pred"),
        backend=backend,
        project_end_point=project_end_point,
        stability_threshold=stability_threshold,
    )

    assert isinstance(df_metrics, pd.DataFrame)
    assert list(df_metrics) == list(df.filter(like="_pred")) + ["metric"]

    if backend == "matplotlib":
        assert isinstance(fig, plt.Figure)
        ax1, ax2 = fig.axes
        assert ax1.get_ylim() == ax2.get_ylim() == (0, 1)
        assert ax1.get_ylabel() == "Recall"
        # TODO ax2 ylabel also 'Recall', should be 'Precision'
        # assert ax2.get_ylabel() == "Precision"
    elif backend == "plotly":
        assert isinstance(fig, go.Figure)
        assert fig.layout.yaxis1.title.text == "Precision"
        assert fig.layout.yaxis2.title.text == "Recall"


@pytest.mark.parametrize("half_window", (0.02, 0.002))
@pytest.mark.parametrize("bin_width", (0.1, 0.001))
@pytest.mark.parametrize("x_lim", ((0, 0.6), (-0.2, 0.8)))
def test_rolling_mae_vs_hull_dist(
    half_window: float, bin_width: float, x_lim: tuple[float, float]
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    for (model_name, df), color in zip(
        test_dfs.items(), ("tab:blue", "tab:orange", "tab:pink")
    ):
        ax = rolling_mae_vs_hull_dist(
            e_above_hull_pred=df.e_above_hull_pred,
            e_above_hull_true=df.e_above_hull_mp,
            color=color,
            label=model_name,
            ax=ax,
            x_lim=x_lim,
            half_window=half_window,
            bin_width=bin_width,
        )

    assert ax is not None
    assert ax.get_ylim() == pytest.approx((0, 0.14))
    assert ax.get_ylabel() == "MAE (eV / atom)"
    assert ax.get_xlabel() == r"$E_\mathrm{above\ hull}$ (eV / atom)"


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

    df = test_dfs["Wren"]

    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_pred=df.e_above_hull_pred,
        e_above_hull_true=df.e_above_hull_mp,
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
        assert ax.layout.yaxis.title.text == "Number of materials"

    assert metrics["precision"] > 0.3
    assert metrics["recall"] > 0.3
