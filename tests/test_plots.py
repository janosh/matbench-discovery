from __future__ import annotations

from typing import Literal

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
    cumulative_clf_metric,
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
    "project_end_point,stability_threshold,metric,expected_line_count",
    [
        ("", 0, "precision", 3),
        ("x", 0, "precision", 6),
        ("x", -0.05, "recall", 6),
        ("xy", 0.1, "recall", 9),
    ],
)
def test_cumulative_precision(
    project_end_point: AxLine,
    stability_threshold: float,
    metric: Literal["precision", "recall"],
    expected_line_count: int,
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    for (model_name, df), color in zip(
        test_dfs.items(), ("tab:blue", "tab:orange", "tab:pink")
    ):
        ax = cumulative_clf_metric(
            e_above_hull_true=df.e_above_hull_mp,
            e_above_hull_pred=df.e_above_hull_pred,
            metric=metric,
            color=color,
            label=model_name,
            project_end_point=project_end_point,
            stability_threshold=stability_threshold,
            ax=ax,
        )

    assert ax is not None
    assert len(ax.lines) == expected_line_count
    assert ax.get_ylim() == (0, 100)
    # assert ax.get_xlim() == pytest.approx((-1.4, 29.4))

    assert ax.get_ylabel() == f"{metric.title()} (%)"


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
