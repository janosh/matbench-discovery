# %%
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from pymatviz.utils import save_fig

from matbench_discovery import PDF_FIGS, SITE_FIGS, plots
from matbench_discovery.preds import df_each_err, models

__author__ = "Janosh Riebesell"
__date__ = "2023-05-25"


# %%
ax = df_each_err[models].plot.box(
    showfliers=False,
    rot=90,
    figsize=(12, 6),
    # color="blue",
    # different fill colors for each box
    # patch_artist=True,
    # notch=True,
    # bootstrap=10000,
    showmeans=True,
    # meanline=True,
)
ax.axhline(0, linewidth=1, color="gray", linestyle="--")


# %%
ax = sns.violinplot(
    data=df_each_err[models], inner="quartile", linewidth=0.3, palette="Set2", width=1
)
ax.set(ylim=(-0.9, 0.9))


# %%
px.violin(
    df_each_err[models].melt(),
    x="variable",
    y="value",
    color="variable",
    violinmode="overlay",
    box=True,
    # points="all",
    hover_data={"variable": False},
    width=1000,
    height=500,
)


# %%
fig = go.Figure()
fig.layout.yaxis.title = plots.quantity_labels["e_above_hull_error"]
fig.layout.margin = dict(l=0, r=0, b=0, t=0)

for col in models:
    val_min = df_each_err[col].quantile(0.05)
    lower_box = df_each_err[col].quantile(0.25)
    median = df_each_err[col].median()
    upper_box = df_each_err[col].quantile(0.75)
    val_max = df_each_err[col].quantile(0.95)

    box_plot = go.Box(
        y=[val_min, lower_box, median, upper_box, val_max],
        name=col,
        width=0.7,
    )
    fig.add_trace(box_plot)

fig.layout.legend.update(orientation="h", y=1.15)
fig.show()


# %%
save_fig(fig, f"{SITE_FIGS}/box-hull-dist-errors.svelte")
save_fig(fig, f"{PDF_FIGS}/box-hull-dist-errors.pdf")
