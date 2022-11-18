# %%
from datetime import datetime

import pandas as pd
import pymatviz

from matbench_discovery import ROOT
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import (
    StabilityCriterion,
    WhichEnergy,
    hist_classified_stable_vs_hull_dist,
    plt,
)

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-08-25"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram is broken down into true positives, false negatives, false positives, and true
negatives based on whether the model predicts candidates to be below the known convex
hull. Ideally, in discovery setting a model should exhibit high recall, i.e. the
majority of materials below the convex hull being correctly identified by the model.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""

today = f"{datetime.now():%Y-%m-%d}"


# %%
dfs = {}
dfs["wren"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")
dfs["m3gnet"] = pd.read_json(
    f"{ROOT}/models/m3gnet/2022-10-31-m3gnet-wbm-IS2RE.json.gz"
).set_index("material_id")
dfs["wrenformer"] = pd.read_csv(
    f"{ROOT}/models/wrenformer/mp/2022-09-20-wrenformer-e_form-ensemble-1-preds.csv"
).set_index("material_id")
dfs["bowsr_megnet"] = pd.read_json(
    f"{ROOT}/models/bowsr/2022-09-22-bowsr-megnet-wbm-IS2RE.json.gz"
).set_index("material_id")


# %%
if "wren" in dfs:
    df = dfs["wren"]
    pred_cols = df.filter(regex=r"_pred_\d").columns
    # make sure we average the expected number of ensemble member predictions
    assert len(pred_cols) == 10
    df["e_form_per_atom_pred"] = df[pred_cols].mean(axis=1)
if "m3gnet" in dfs:
    df = dfs["m3gnet"]
    df["e_form_per_atom_pred"] = df.e_form_per_atom_m3gnet
if "bowsr_megnet" in dfs:
    df = dfs["bowsr_megnet"]
    df["e_form_per_atom_pred"] = df.e_form_per_atom_bowsr


# %%
which_energy: WhichEnergy = "true"
stability_crit: StabilityCriterion = "energy"
fig, axs = plt.subplots(2, 3, figsize=(18, 9))

model_name = "m3gnet"
df = dfs[model_name]

df["e_above_hull_mp"] = df_wbm.e_above_hull_mp2020_corrected_ppd_mp
df["e_form_per_atom"] = df_wbm.e_form_per_atom_mp2020_corrected


for batch_idx, ax in zip(range(1, 6), axs.flat):
    batch_df = df[df.index.str.startswith(f"wbm-step-{batch_idx}-")]
    assert 1e4 < len(batch_df) < 1e5, print(f"{len(batch_df) = :,}")

    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_pred=batch_df.e_form_per_atom_pred - batch_df.e_form_per_atom,
        e_above_hull_true=batch_df.e_above_hull_mp,
        which_energy=which_energy,
        stability_crit=stability_crit,
        ax=ax,
    )

    text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"
    ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

    title = f"Batch {batch_idx} ({len(batch_df.filter(like='e_').dropna()):,})"
    ax.set(title=title)


ax, metrics = hist_classified_stable_vs_hull_dist(
    e_above_hull_pred=df.e_form_per_atom_pred - df.e_form_per_atom,
    e_above_hull_true=df.e_above_hull_mp,
    which_energy=which_energy,
    stability_crit=stability_crit,
    ax=axs.flat[-1],
)

text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"
ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

axs.flat[-1].set(title=f"Combined ({len(df.filter(like='e_').dropna()):,})")
axs.flat[0].legend(frameon=False, loc="upper left")

img_name = f"{today}-{model_name}-wbm-hull-dist-hist-batches"
suptitle = img_name.replace("-", "/", 2).replace("-", " ")
fig.suptitle(suptitle, y=1.07, fontsize=16)


# %%
ax.figure.savefig(f"{ROOT}/figures/{img_name}.pdf")


# %%
pymatviz.density_scatter(
    df=dfs["m3gnet"].query("e_form_per_atom < 5"),
    x="e_form_per_atom",
    y="e_form_per_atom_pred",
)
