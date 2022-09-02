# %%
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import precision_recall_vs_calc_count


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=200)
plt.rc("font", size=16)


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

dfs: dict[str, pd.DataFrame] = {}
for model_name in ("Wren", "CGCNN", "Voronoi"):
    dfs[model_name] = pd.read_csv(
        f"{ROOT}/data/2022-06-11-from-rhys/{model_name}-mp-initial-structures.csv"
    ).set_index("material_id")

# dfs["M3GNet"] = pd.read_json(
#     f"{ROOT}/data/2022-08-16-m3gnet-wbm-relax-results-IS2RE.json.gz"
# ).set_index("material_id")

dfs["Wrenformer"] = pd.read_csv(
    f"{ROOT}/data/2022-08-16-wrenformer-ensemble-predictions.csv.bz2"
).set_index("material_id")


# download wbm-steps-summary.csv (23.31 MB)
df_summary = pd.read_csv(
    "https://figshare.com/ndownloader/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")


# %%
for (model_name, df), color in zip(
    dfs.items(),
    ("tab:blue", "tab:orange", "teal", "tab:pink", "black", "red", "turquoise"),
):
    df["e_above_mp_hull"] = df_hull.e_above_mp_hull

    assert df.e_above_mp_hull.isna().sum() == 0

    target_col = "e_form_target"
    rare = "all"

    # from pymatgen.core import Composition
    # rare = "no-lanthanides"
    # df["contains_rare_earths"] = df.composition.map(
    #     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
    # )
    # df = df.query("~contains_rare_earths")

    e_above_mp_hull = df.e_above_mp_hull

    try:
        if model_name == "M3GNet":
            model_preds = df.e_form_m3gnet
            targets = df.e_form_wbm
        elif "Wrenformer" in model_name:
            df["e_form_per_atom_pred_ens"] = df.e_form_pred_ens / df.n_sites
            df["e_form_per_atom"] = df.e_form / df.n_sites
            model_preds = df.e_form_per_atom_pred_ens
            targets = df.e_form_per_atom
        elif df.filter(regex=r"_pred_\d").shape[1] > 1:
            assert df.filter(regex=r"_pred_\d").shape[1] == 10
            model_preds = df.filter(regex=r"_pred_\d").mean(axis=1)
            targets = df.e_form_target
        elif "e_form_pred" in df and "e_form_target" in df:
            model_preds = df.e_form_pred
            targets = df.e_form_target
        else:
            raise ValueError(f"Unhandled {model_name = }")
    except AttributeError as exc:
        raise KeyError(f"{model_name = }") from exc

    df["residual"] = model_preds - targets + df.e_above_mp_hull

    ax = precision_recall_vs_calc_count(
        df,
        residual_col="residual",
        e_above_hull_col="e_above_mp_hull",
        color=color,
        label=model_name,
        intersect_lines="recall_xy",
        # intersect_lines="all",
    )

ax.legend(frameon=False, loc="lower right")

ax.figure.set_size_inches(10, 9)

img_path = f"{ROOT}/figures/{today}-precision-recall-vs-calc-count-{rare=}.pdf"
# plt.savefig(img_path)
