"""Compare CHGNet long vs short relaxations."""

# %%
import numpy as np
import pandas as pd
import pymatviz as pmv
from aviary.wren.utils import get_prototype_from_protostructure
from IPython.display import display
from pymatviz.enums import Key
from pymatviz.io import df_to_html, df_to_pdf, save_fig
from pymatviz.powerups import add_identity_line
from pymatviz.ptable import ptable_heatmap_plotly
from pymatviz.utils import PLOTLY, bin_df_cols

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.data import DataFiles, Model, df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.preds import df_each_pred, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-03-20"


# %%
model = Model.wrenformer.label
model_low = model.lower()
max_each_true = 1
min_each_pred = 1
df_each_pred[MbdKey.each_true] = df_preds[MbdKey.each_true]
bad_ids = df_each_pred.query(
    f"{model} > {min_each_pred} & {MbdKey.each_true} < {min_each_pred}"
).index

df_wbm[Key.spg_num] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_bad = df_wbm.loc[bad_ids]
title = f"{len(df_bad)} {model} preds<br>with {max_each_true=}, {min_each_pred=}"


# %%
df_mp = pd.read_csv(DataFiles.mp_energies.path).set_index(Key.mat_id)
df_mp[Key.spg_num] = df_mp[Key.wyckoff].str.split("_").str[2].astype(int)
df_mp["isopointal_proto_from_aflow"] = df_mp[Key.wyckoff].map(
    get_prototype_from_protostructure
)
df_mp.isopointal_proto_from_aflow.value_counts().head(12)


# %%
fig = pmv.spacegroup_bar(df_bad[Key.spg_num])
fig.layout.title.update(text=f"Spacegroup hist for {title}", y=0.96)
fig.layout.margin.update(l=0, r=0, t=80, b=0)
save_fig(fig, f"{PDF_FIGS}/spacegroup-hist-{model.lower()}-failures.pdf")
fig.show()


# %%
proto_col = "Isopointal Prototypes"
df_proto_counts = (
    df_bad[MbdKey.init_wyckoff]
    .map(get_prototype_from_protostructure)
    .value_counts()
    .to_frame()
)


df_proto_counts["MP occurrences"] = 0
mp_proto_counts = df_mp.isopointal_proto_from_aflow.value_counts()
for proto in df_proto_counts.index:
    df_proto_counts.loc[proto, "MP occurrences"] = mp_proto_counts.get(proto, 0)

df_proto_counts = df_proto_counts.reset_index(names=proto_col)

# improve proto_col readability
df_proto_counts[proto_col] = df_proto_counts[proto_col].str.replace("_", "-")

styler = df_proto_counts.head(10).style.background_gradient(cmap="viridis")
styler.set_caption(f"Top 10 {proto_col} in {len(df_bad)} {model} failures")
display(styler)
img_name = f"proto-counts-{model_low}-failures"
df_to_html(styler, file_path=f"{SITE_FIGS}/{img_name}.svelte")
df_to_pdf(styler, f"{PDF_FIGS}/{img_name}.pdf")


# %%
fig = pmv.spacegroup_sunburst(
    df_bad[Key.spg_num], width=350, height=350, show_counts="percent"
)
# fig.layout.title.update(text=f"Spacegroup sunburst for {title}", x=0.5, font_size=14)
fig.layout.margin.update(l=1, r=1, t=1, b=1)
fig.show()


# %%
save_fig(fig, f"{PDF_FIGS}/spacegroup-sunburst-{model_low}-failures.pdf")
save_fig(fig, f"{SITE_FIGS}/spacegroup-sunburst-{model_low}-failures.svelte")


# %%
fig = ptable_heatmap_plotly(df_bad[Key.formula])
fig.layout.title = f"Elements in {title}"
fig.layout.margin = dict(l=0, r=0, t=50, b=0)
fig.show()
save_fig(fig, f"{PDF_FIGS}/elements-{model_low}-failures.pdf")


# %%
model = Model.wrenformer.label
cols = [model, MbdKey.each_true]
bin_cnt_col = "bin counts"
df_bin = bin_df_cols(
    df_each_pred, [MbdKey.each_true, model], n_bins=200, bin_counts_col=bin_cnt_col
)
log_cnt_col = f"log {bin_cnt_col}"
df_bin[log_cnt_col] = np.log1p(df_bin[bin_cnt_col]).round(2)


# %%
fig = df_bin.reset_index().plot.scatter(
    x=MbdKey.each_true,
    y=model,
    hover_data=cols,
    hover_name=df_preds.index.name,
    backend=PLOTLY,
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
img_name = "hull-dist-parity-wrenformer-failures"
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=600, height=300)
