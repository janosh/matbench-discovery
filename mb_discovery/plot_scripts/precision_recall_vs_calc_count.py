# %%
from datetime import datetime

import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts import plt
from mb_discovery.plot_scripts.plot_funcs import (
    StabilityCriterion,
    precision_recall_vs_calc_count,
)


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"


# %%
DATA_DIR = f"{ROOT}/data/2022-06-11-from-rhys"
df_hull = pd.read_csv(f"{DATA_DIR}/wbm-e-above-mp-hull.csv").set_index("material_id")

dfs: dict[str, pd.DataFrame] = {}
for model_name in ("Wren", "CGCNN", "Voronoi"):
    df = pd.read_csv(
        f"{DATA_DIR}/{model_name.lower()}-mp-initial-structures.csv"
    ).set_index("material_id")
    dfs[model_name] = df

# dfs["M3GNet"] = pd.read_json(
#     f"{ROOT}/data/2022-08-16-m3gnet-wbm-relax-results-IS2RE.json.gz"
# ).set_index("material_id")

dfs["Wrenformer"] = pd.read_csv(
    f"{ROOT}/data/2022-08-16-wrenformer-preds.csv.bz2"
).set_index("material_id")


# %% download wbm-steps-summary.csv (23.31 MB)
df_summary = pd.read_csv(
    "https://figshare.com/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")


stability_crit: StabilityCriterion = "energy"


# %%
for (model_name, df), color in zip(
    dfs.items(),
    ("tab:blue", "tab:orange", "teal", "tab:pink", "black", "red", "turquoise"),
):
    rare = "all"

    # from pymatgen.core import Composition
    # rare = "no-lanthanides"
    # df["contains_rare_earths"] = df.composition.map(
    #     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
    # )
    # df = df.query("~contains_rare_earths")

    if "std" in stability_crit:
        # TODO column names to compute standard deviation from are currently hardcoded
        # needs to be updated when adding non-aviary models with uncertainty estimation
        var_aleatoric = (df.filter(regex=r"_ale_\d") ** 2).mean(axis=1)
        var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
        std_total = (var_epistemic + var_aleatoric) ** 0.5
    else:
        std_total = None

    try:
        if model_name == "M3GNet":
            model_preds = df.e_form_m3gnet
            targets = df.e_form_wbm
        elif "Wrenformer" in model_name:
            df["e_form_per_atom_pred_ens"] = df.e_form_pred_ens / df.n_sites
            df["e_form_per_atom"] = df.e_form / df.n_sites
            model_preds = df.e_form_per_atom_pred_ens
            targets = df.e_form_per_atom
        elif len(pred_cols := df.filter(like="e_form_pred").columns) >= 1:
            # Voronoi+RF has single prediction column, Wren and CGCNN each have 10
            # other cases are unexpected
            assert len(pred_cols) in (1, 10), f"{model_name=} has {len(pred_cols)=}"
            model_preds = df[pred_cols].mean(axis=1)
            targets = df.e_form_target
        else:
            raise ValueError(f"Unhandled {model_name = }")
    except AttributeError as exc:
        raise KeyError(f"{model_name = }") from exc

    df["e_above_mp_hull"] = df_hull.e_above_mp_hull
    df["e_above_hull_pred"] = model_preds - targets

    ax = precision_recall_vs_calc_count(
        e_above_hull_error=df.e_above_hull_pred + df.e_above_mp_hull,
        e_above_hull_true=df.e_above_mp_hull,
        color=color,
        label=model_name,
        intersect_lines="recall_xy",  # or "precision_xy", None, 'all'
        stability_crit=stability_crit,
        std_pred=std_total,
    )

ax.legend(frameon=False, loc="lower right")

ax.figure.set_size_inches(10, 9)

img_path = f"{ROOT}/figures/{today}-precision-recall-vs-calc-count-{rare=}.pdf"
plt.savefig(img_path)
