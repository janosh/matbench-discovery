from __future__ import annotations

from typing import Any, Literal

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from matbench_discovery import ROOT
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import (
    AxLine,
    StabilityCriterion,
    cumulative_clf_metric,
    hist_classified_stable_as_func_of_hull_dist,
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

    df["e_above_hull_pred"] = model_preds - df.e_form_target

    test_dfs[model_name] = df


@pytest.mark.parametrize(
    "project_end_point,stability_crit,stability_threshold,metric,expected_line_count",
    [
        ("", "energy", 0, "precision", 3),
        ("x", "energy-std", 0, "precision", 6),
        ("x", "energy", -0.05, "recall", 6),
        ("xy", "energy-std", 0.1, "recall", 9),
        # TODO: test energy+std
    ],
)
def test_cumulative_precision(
    project_end_point: AxLine,
    stability_crit: StabilityCriterion,
    stability_threshold: float,
    metric: Literal["precision", "recall"],
    expected_line_count: int,
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    for (model_name, df), color in zip(
        test_dfs.items(), ("tab:blue", "tab:orange", "tab:pink")
    ):
        if "std" in stability_crit:
            std_pred = df.filter(like=r"_pred").std(axis=1)
        else:
            std_pred = None

        ax = cumulative_clf_metric(
            e_above_hull_error=df.e_above_hull_pred,
            e_above_hull_true=df.e_above_hull_mp,
            metric=metric,
            color=color,
            label=model_name,
            project_end_point=project_end_point,
            stability_crit=stability_crit,
            std_pred=std_pred,
            stability_threshold=stability_threshold,
            ax=ax,
        )

    assert ax is not None
    assert len(ax.lines) == expected_line_count
    assert ax.get_ylim() == (0, 100)
    # assert ax.get_xlim() == pytest.approx((-1.4, 29.4))

    assert ax.get_ylabel() == f"{metric.title()} (%)"


@pytest.mark.parametrize(
    "kwargs, expected_exc, match_pat",
    [
        (dict(stability_crit="INVALID"), ValueError, "Invalid stability_crit="),
    ],
)
def test_cumulative_precision_raises(
    kwargs: dict[str, Any], expected_exc: type[Exception], match_pat: str
) -> None:
    with pytest.raises(expected_exc, match=match_pat):
        cumulative_clf_metric(
            e_above_hull_error=test_dfs["Wren"].e_above_hull_pred,
            e_above_hull_true=test_dfs["Wren"].e_above_hull_mp,
            metric="precision",
            **kwargs,
        )


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
    assert ax.get_xlabel() == r"$\Delta E_{Hull-MP}$ (eV / atom)"


@pytest.mark.parametrize("stability_threshold", (0.1, 0.01))
@pytest.mark.parametrize("stability_crit", ("energy", "energy+std", "energy-std"))
@pytest.mark.parametrize("x_lim", ((0, 0.6), (-0.2, 0.8)))
def test_hist_classified_stable_as_func_of_hull_dist(
    stability_threshold: float,
    stability_crit: StabilityCriterion,
    x_lim: tuple[float, float],
) -> None:
    ax = plt.figure().gca()  # new figure ensures test functions use different axes

    df = test_dfs["Wren"]

    if "std" in stability_crit:
        var_aleatoric = (df.filter(like="_ale_") ** 2).mean(axis=1)
        var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
        std_total = (var_epistemic + var_aleatoric) ** 0.5
    else:
        std_total = None

    ax, metrics = hist_classified_stable_as_func_of_hull_dist(
        e_above_hull_pred=df.e_above_hull_pred,
        e_above_hull_true=df.e_above_hull_mp,
        ax=ax,
        stability_threshold=stability_threshold,
        stability_crit=stability_crit,
        x_lim=x_lim,
        std_pred=std_total,
    )

    assert ax is not None
    # assert ax.get_ylim() == pytest.approx((0, 6.3))
    assert ax.get_ylabel() == "Number of compounds"
    assert ax.get_xlabel() == r"$\Delta E_{Hull-MP}$ (eV / atom)"

    assert metrics["precision"] > 0.3
    assert metrics["recall"] > 0.3
