"""Compare CHGNet long vs short relaxations."""


# %%
from pymatviz import spacegroup_hist, spacegroup_sunburst
from pymatviz.ptable import ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import ROOT
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import df_each_pred, df_preds, each_true_col

__author__ = "Janosh Riebesell"
__date__ = "2023-03-20"


# %%
model = "Wrenformer"
max_each_true = 1
min_each_pred = 1
df_each_pred[each_true_col] = df_preds[each_true_col]
bad_ids = df_each_pred.query(
    f"{model} > {min_each_pred} & {each_true_col} < {min_each_pred}"
).index

spg_col = "spacegroup"
wyk_col = "wyckoff_spglib"
df_wbm[spg_col] = df_wbm[wyk_col].str.split("_").str[2].astype(int)
df_bad = df_wbm.loc[bad_ids]
title = f"{len(df_bad)} {model} preds with {max_each_true=}, {min_each_pred=}"


# %%
ax = spacegroup_hist(df_bad[spg_col])
fig = spacegroup_sunburst(df_bad[spg_col])
fig.layout.title = f"Spacegroup sunburst for {title}"
ax.set_title(f"Spacegroup hist for {title}", y=1.15)
out_dir = f"{ROOT}/tmp/figures"
save_fig(fig, f"{out_dir}/bad-{model}-spacegroup-sunburst.png", scale=3)
save_fig(ax, f"{out_dir}/bad-{model}-spacegroup-hist.png", dpi=300)


# %%
fig = ptable_heatmap_plotly(df_bad.formula)
fig.layout.title = f"Elements in {title}"
fig.layout.margin = dict(l=0, r=0, t=50, b=0)
save_fig(fig, f"{out_dir}/bad-{model}-elements.png", scale=3)
