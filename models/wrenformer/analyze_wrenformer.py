"""Compare CHGNet long vs short relaxations."""


# %%
import numpy as np
import pandas as pd
from aviary.wren.utils import get_isopointal_proto_from_aflow
from pymatviz import spacegroup_hist, spacegroup_sunburst
from pymatviz.io import df_to_pdf, df_to_svelte_table
from pymatviz.ptable import ptable_heatmap_plotly
from pymatviz.utils import add_identity_line, bin_df_cols, save_fig

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.data import DATA_FILES, df_wbm
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
wyckoff_col = "wyckoff_spglib"
df_wbm[spg_col] = df_wbm[wyckoff_col].str.split("_").str[2].astype(int)
df_bad = df_wbm.loc[bad_ids]
title = f"{len(df_bad)} {model} preds<br>with {max_each_true=}, {min_each_pred=}"


# %%
df_mp = pd.read_csv(DATA_FILES.mp_energies).set_index("material_id")
df_mp[spg_col] = df_mp[wyckoff_col].str.split("_").str[2].astype(int)
df_mp["isopointal_proto_from_aflow"] = df_mp[wyckoff_col].map(
    get_isopointal_proto_from_aflow
)
df_mp.isopointal_proto_from_aflow.value_counts().head(12)


# %%
ax = spacegroup_hist(df_bad[spg_col])
ax.set_title(f"Spacegroup hist for {title}", y=1.15)
save_fig(ax, f"{PDF_FIGS}/spacegroup-hist-{model.lower()}-failures.pdf")


# %%
proto_col = "Isopointal Prototypes"
df_proto_counts = (
    df_bad[wyckoff_col].map(get_isopointal_proto_from_aflow).value_counts().to_frame()
)


df_proto_counts["MP occurrences"] = 0
mp_proto_counts = df_mp.isopointal_proto_from_aflow.value_counts()
for proto in df_proto_counts.index:
    df_proto_counts.loc[proto, "MP occurrences"] = mp_proto_counts.get(proto, 0)

df_proto_counts = df_proto_counts.reset_index(names=proto_col)

# improve proto_col readability
df_proto_counts[proto_col] = df_proto_counts[proto_col].str.replace("_", "-")

styler = df_proto_counts.head(10).style.background_gradient(cmap="viridis")
styles = {
    "": "font-family: sans-serif; border-collapse: collapse;",
    "td, th": "border: none; padding: 4px 6px; white-space: nowrap;",
    "th": "border: 1px solid; border-width: 1px 0; text-align: left;",
}
styler.set_table_styles([dict(selector=sel, props=styles[sel]) for sel in styles])
styler.set_uuid("")

df_to_svelte_table(styler, f"{SITE_FIGS}/proto-counts-{model}-failures.svelte")
df_to_pdf(styler, f"{PDF_FIGS}/proto-counts-{model}-failures.pdf")


# %%
fig = spacegroup_sunburst(df_bad[spg_col], width=350, height=350)
# fig.layout.title.update(text=f"Spacegroup sunburst for {title}", x=0.5, font_size=14)
fig.layout.margin.update(l=1, r=1, t=1, b=1)
fig.show()


# %%
save_fig(fig, f"{PDF_FIGS}/spacegroup-sunburst-{model.lower()}-failures.pdf")
save_fig(fig, f"{SITE_FIGS}/spacegroup-sunburst-{model}-failures.svelte")


# %%
fig = ptable_heatmap_plotly(df_bad.formula)
fig.layout.title = f"Elements in {title}"
fig.layout.margin = dict(l=0, r=0, t=50, b=0)
fig.show()
save_fig(fig, f"{PDF_FIGS}/elements-{model.lower()}-failures.pdf")


# %%
model = "Wrenformer"
cols = [model, each_true_col]
bin_cnt_col = "bin counts"
df_bin = bin_df_cols(
    df_each_pred, [each_true_col, model], n_bins=200, bin_counts_col=bin_cnt_col
)
log_cnt_col = f"log {bin_cnt_col}"
df_bin[log_cnt_col] = np.log1p(df_bin[bin_cnt_col]).round(2)


# %%
fig = df_bin.reset_index().plot.scatter(
    x=each_true_col,
    y=model,
    hover_data=cols,
    hover_name=df_preds.index.name,
    backend="plotly",
    color=log_cnt_col,
    color_continuous_scale="turbo",
)

# title = "Analysis of Wrenformer failure cases in the highlighted rectangle"
# fig.layout.title.update(text=title, x=0.5)
fig.layout.margin.update(l=0, r=0, t=0, b=0)
fig.layout.legend.update(title="", x=1, y=0, xanchor="right")
add_identity_line(fig)
fig.layout.coloraxis.colorbar.update(
    x=1, y=0.5, xanchor="right", thickness=12, title=""
)
# add shape shaded rectangle at x < 1, y > 1
fig.add_shape(
    type="rect", **dict(x0=1, y0=1, x1=-1, y1=6), fillcolor="gray", opacity=0.2
)
fig.show()


# %%
img_name = "hull-dist-scatter-wrenformer-failures"
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=600, height=300)
