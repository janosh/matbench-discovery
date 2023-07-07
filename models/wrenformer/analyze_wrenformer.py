"""Compare CHGNet long vs short relaxations."""


# %%
import pandas as pd
from aviary.wren.utils import get_isopointal_proto_from_aflow
from pymatviz import spacegroup_hist, spacegroup_sunburst
from pymatviz.ptable import ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, PDF_FIGS
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.plots import df_to_pdf, df_to_svelte_table
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

df_to_svelte_table(styler, f"{FIGS}/proto-counts-{model}-failures.svelte")
df_to_pdf(styler, f"{PDF_FIGS}/proto-counts-{model}-failures.pdf")


# %%
fig = spacegroup_sunburst(df_bad[spg_col], width=350, height=350)
fig.layout.title.update(text=f"Spacegroup sunburst for {title}", x=0.5, font_size=14)
fig.show()
save_fig(fig, f"{PDF_FIGS}/spacegroup-sunburst-{model.lower()}-failures.pdf")
# save_fig(fig, f"{FIGS}/spacegroup-sunburst-{model}-failures.svelte")


# %%
fig = ptable_heatmap_plotly(df_bad.formula)
fig.layout.title = f"Elements in {title}"
fig.layout.margin = dict(l=0, r=0, t=50, b=0)
fig.show()
save_fig(fig, f"{PDF_FIGS}/elements-{model.lower()}-failures.pdf")
