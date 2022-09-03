from __future__ import annotations

from typing import Any, Sequence

import pandas as pd
import pytest

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import precision_recall_vs_calc_count


df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

test_dfs: dict[str, pd.DataFrame] = {}
for model_name in ("Wren", "CGCNN", "Voronoi"):
    df = pd.read_csv(
        f"{ROOT}/data/2022-06-11-from-rhys/{model_name}-mp-initial-structures.csv",
        nrows=100,
    ).set_index("material_id")

    df["e_above_mp_hull"] = df_hull.e_above_mp_hull

    test_dfs[model_name] = df


@pytest.mark.parametrize(
    "intersect_lines, stability_crit, stability_threshold, expected_line_count",
    [
        ((), "energy", 0, 11),
        ("precision_x", "energy+std", 0, 23),
        (["recall_y"], "energy", -0.1, 35),
        ("all", "energy-std", 0.1, 56),
    ],
)
def test_precision_recall_vs_calc_count(
    intersect_lines: str | Sequence[str],
    stability_crit: str,
    stability_threshold: float,
    expected_line_count: int,
) -> None:
    ax = None

    for (model_name, df), color in zip(
        test_dfs.items(), ("tab:blue", "tab:orange", "tab:pink")
    ):
        model_preds = df.filter(like=r"_pred").mean(axis=1)
        targets = df.e_form_target

        df["residual"] = model_preds - targets + df.e_above_mp_hull

        ax = precision_recall_vs_calc_count(
            df,
            residual_col="residual",
            e_above_hull_col="e_above_mp_hull",
            color=color,
            label=model_name,
            intersect_lines=intersect_lines,
            stability_crit=stability_crit,  # type: ignore[arg-type]
            stability_threshold=stability_threshold,
            ax=ax,
        )

    assert ax is not None
    assert len(ax.lines) == expected_line_count
    assert ax.get_ylim() == (0, 100)
    assert ax.get_xlim() == pytest.approx((-1.4, 29.4))


@pytest.mark.parametrize(
    "kwargs, expected_exc, match_pat",
    [
        (dict(intersect_lines="INVALID"), ValueError, "Invalid intersect_lines="),
        (dict(stability_crit="INVALID"), ValueError, "Invalid stability_crit="),
    ],
)
def test_precision_recall_vs_calc_count_raises(
    kwargs: dict[str, Any], expected_exc: type[Exception], match_pat: str
) -> None:
    with pytest.raises(expected_exc, match=match_pat):
        precision_recall_vs_calc_count(
            test_dfs["Wren"],
            residual_col="residual",
            e_above_hull_col="e_above_mp_hull",
            **kwargs,
        )
