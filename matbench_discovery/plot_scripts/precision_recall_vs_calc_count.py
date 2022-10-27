# %%
from datetime import datetime

import pandas as pd
from sklearn.metrics import f1_score

from matbench_discovery import ROOT
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import StabilityCriterion, precision_recall_vs_calc_count

__author__ = "Rhys Goodall, Janosh Riebesell"

today = f"{datetime.now():%Y-%m-%d}"


# %%
rare = "all"

dfs: dict[str, pd.DataFrame] = {}
for model_name in ("wren", "cgcnn", "voronoi"):
    csv_path = (
        f"{ROOT}/data/2022-06-11-from-rhys/{model_name}-mp-initial-structures.csv"
    )
    df = pd.read_csv(csv_path).set_index("material_id")
    dfs[model_name] = df

dfs["m3gnet"] = pd.read_json(
    f"{ROOT}/models/m3gnet/2022-08-16-m3gnet-wbm-IS2RE.json.gz"
).set_index("material_id")

dfs["wrenformer"] = pd.read_csv(
    f"{ROOT}/models/wrenformer/mp/"
    "2022-09-20-wrenformer-e_form-ensemble-1-preds-e_form_per_atom.csv"
).set_index("material_id")

dfs["bowsr_megnet"] = pd.read_json(
    f"{ROOT}/models/bowsr/2022-09-22-bowsr-megnet-wbm-IS2RE.json.gz"
).set_index("material_id")

print(f"loaded models: {list(dfs)}")


# %%
stability_crit: StabilityCriterion = "energy"
colors = "tab:blue tab:orange teal tab:pink black red turquoise tab:purple".split()
F1s: dict[str, float] = {}

for model_name, df in dfs.items():
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
        if model_name == "m3gnet":
            model_preds = df.e_form_m3gnet
        elif "wrenformer" in model_name:
            model_preds = df.e_form_per_atom_pred_ens
        elif len(pred_cols := df.filter(like="e_form_pred").columns) >= 1:
            # Voronoi+RF has single prediction column, Wren and CGCNN each have 10
            # other cases are unexpected
            assert len(pred_cols) in (1, 10), f"{model_name=} has {len(pred_cols)=}"
            model_preds = df[pred_cols].mean(axis=1)
        elif "bowsr" in model_name:
            model_preds = df.e_form_per_atom_bowsr
        else:
            raise ValueError(f"Unhandled {model_name = }")
    except AttributeError as exc:
        raise KeyError(f"{model_name = }") from exc

    df["e_above_hull_mp"] = df_wbm.e_above_hull_mp2020_corrected
    df["e_form_per_atom"] = df_wbm.e_form_per_atom_mp2020_corrected
    df["e_above_hull_pred"] = model_preds - df.e_form_per_atom
    if n_nans := df.isna().values.sum() > 0:
        assert n_nans < 10, f"{model_name=} has {n_nans=}"
        df = df.dropna()

    F1 = f1_score(df.e_above_hull_mp < 0, df.e_above_hull_pred < 0)
    F1s[model_name] = F1


# %%
for (model_name, F1), color in zip(sorted(F1s.items(), key=lambda x: x[1]), colors):
    df = dfs[model_name]

    ax = precision_recall_vs_calc_count(
        e_above_hull_error=df.e_above_hull_pred + df.e_above_hull_mp,
        e_above_hull_true=df.e_above_hull_mp,
        color=color,
        label=f"{model_name} {F1=:.2}",
        intersect_lines="recall_xy",  # or "precision_xy", None, 'all'
        stability_crit=stability_crit,
        std_pred=std_total,
    )

# optimal recall line finds all stable materials without any false positives
# can be included to confirm all models start out of with near optimal recall
# and to see how much each model overshoots total n_stable
n_below_hull = sum(df_wbm.e_above_hull_mp2020_corrected < 0)
ax.plot(
    [0, n_below_hull],
    [0, 100],
    color="green",
    linestyle="dashed",
    linewidth=1,
    label="Optimal Recall",
)

ax.figure.set_size_inches(10, 9)
ax.set(xlim=(0, None))
# keep this outside loop so all model names appear in legend
ax.legend(frameon=False, loc="lower right")

img_name = f"{today}-precision-recall-vs-calc-count-{rare=}"
ax.set(title=img_name.replace("-", "/", 2).replace("-", " ").title())
# x-ticks every 10k materials
ax.set(xticks=range(0, int(ax.get_xlim()[1]), 10_000))


# %%
ax.figure.savefig(f"{ROOT}/figures/{img_name}.pdf")
