# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from mb_discovery import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
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

dfs["M3GNet"] = pd.read_json(
    f"{ROOT}/data/2022-08-16-m3gnet-wbm-relax-results-IS2RE.json.gz"
).set_index("material_id")

dfs["Wrenformer"] = pd.read_csv(
    f"{ROOT}/data/2022-08-16-wrenformer-ensemble-predictions.csv.bz2"
).set_index("material_id")

# dfs["Wrenformer"]["e_form_target"] = dfs["Wren"]["e_form_target"]
# dfs["M3GNet"]["e_form_target"] = dfs["Wren"]["e_form_target"]


# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

for model_name, color in zip(
    ("Wren", "CGCNN", "Voronoi", "M3GNet", "Wrenformer"),
    ("tab:blue", "tab:orange", "teal", "tab:pink", "black"),
    strict=True,
):
    df = dfs[model_name]
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
        elif model_name == "Wrenformer":
            model_preds = df.e_form_pred_ens
            targets = df.e_form
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
    df = df.sort_values(by="residual")

    # epistemic_var = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)

    # aleatoric_var = (df.filter(like="_ale_") ** 2).mean(axis=1)

    # std_total = (epistemic_var + aleatoric_var) ** 0.5

    # criterion = "std"
    # test = df.residual + std_total

    # criterion = "neg"
    # test = df.residual - std_total

    criterion = "energy"

    # stability_thresh = 0.02
    stability_thresh = 0
    # stability_thresh = 0.10

    true_pos_mask = (df.e_above_mp_hull <= stability_thresh) & (
        df.residual <= stability_thresh
    )
    false_neg_mask = (df.e_above_mp_hull <= stability_thresh) & (
        df.residual > stability_thresh
    )
    false_pos_mask = (df.e_above_mp_hull > stability_thresh) & (
        df.residual <= stability_thresh
    )

    energy_type = "pred"
    true_pos_cumsum = true_pos_mask.cumsum()
    xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

    ppv = true_pos_cumsum / (true_pos_cumsum + false_pos_mask.cumsum()) * 100
    n_true_pos = sum(true_pos_mask)
    n_false_neg = sum(false_neg_mask)
    n_total_pos = n_true_pos + n_false_neg
    tpr = true_pos_cumsum / n_total_pos * 100

    end = int(np.argmax(tpr))

    xs = np.arange(end)

    precision_curve = interp1d(xs, ppv[:end], kind="cubic")
    rolling_recall_curve = interp1d(xs, tpr[:end], kind="cubic")

    line_kwargs = dict(
        linewidth=3,
        color=color,
        markevery=[-1],
        marker="x",
        markersize=14,
        markeredgewidth=2.5,
    )
    ax.plot(xs, precision_curve(xs), linestyle="-", **line_kwargs)
    ax.plot(xs, rolling_recall_curve(xs), linestyle=":", **line_kwargs)
    ax.plot((0, 0), (0, 0), label=model_name, **line_kwargs)


ax.set(xlabel="Number of Calculations", ylabel="Percentage")

ax.set(xlim=(0, 8e4), ylim=(0, 100))

model_legend = ax.legend(frameon=False, loc="lower right")
ax.add_artist(model_legend)

[precision] = ax.plot((0, 0), (0, 0), "black", linestyle="-")
[recall] = ax.plot((0, 0), (0, 0), "black", linestyle=":")
ax.legend(
    [precision, recall], ("Precision", "Recall"), frameon=False, loc="upper right"
)

img_path = (
    f"{ROOT}/figures/{today}-precision-recall-vs-calc-count-"
    f"{energy_type=}-{criterion=}-{rare=}.pdf"
)
# plt.savefig(img_path)
