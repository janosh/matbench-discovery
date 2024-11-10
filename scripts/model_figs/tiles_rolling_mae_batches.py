"""Plot rolling MAE as a function of hull distance for a single model but split per WBM
batch in a single plot.
"""

# %%
import math

import pymatviz as pmv
from plotly.subplots import make_subplots
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds.discovery import df_each_pred, df_preds
from matbench_discovery.preds.discovery import models as all_models

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

batch_col = "batch_idx"
df_each_pred[batch_col] = "Batch " + df_each_pred.index.str.split("-").str[1]
df_err, df_std = None, None  # variables to cache rolling MAE and std
models = globals().get("models", all_models)


save_individual_figs = globals().get("save_individual_figs", True)
test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]

n_cols = 3
use_full_rows = globals().get("use_full_rows", True)
if use_full_rows:
    # drop last models that don't fit in last row
    n_rows = len(models_to_plot) // n_cols
    models_to_plot = models_to_plot[: n_rows * n_cols]
else:
    n_rows = math.ceil(len(models) / n_cols)


# %%
# Create subplots with one row per column in the DataFrame
fig = make_subplots(
    rows=n_rows,
    cols=n_cols,
    subplot_titles=models_to_plot,
    shared_xaxes=True,
    shared_yaxes=True,
    vertical_spacing=0.04,
    horizontal_spacing=0.03,
)
# Update title font size for all subplot titles
fig.layout.update(height=230 * n_rows)
fig.layout.update(width=280 * n_cols)

subfig = None
for i, model in enumerate(models_to_plot):
    df_pivot = df_each_pred.pivot(columns=batch_col, values=model)  # noqa: PD010

    subfig, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_preds[MbdKey.each_true],
        e_above_hull_preds=df_pivot,
        # df_rolling_err=df_err,
        # df_err_std=df_std,
        show_dummy_mae=False,
        with_sem=False,
    )

    if save_individual_figs:
        model_snake_case = model.lower().replace(" + ", "-").replace(" ", "-")
        img_path = f"rolling-mae-vs-hull-dist-wbm-batches-{model_snake_case}"
        subfig.layout.margin.update(l=10, r=10, b=10, t=10)
        subfig.layout.legend.update(title=f"<b>{model}</b>")
        subfig.layout.update(hovermode="x unified", hoverlabel_bgcolor="black")

        subfig.update_traces(
            hovertemplate="y=%{y:.3f} eV",
            selector=lambda trace: trace.name.startswith("Batch"),
        )
        subfig.show()
        pmv.save_fig(subfig, f"{SITE_FIGS}/{img_path}.svelte")
        pmv.save_fig(subfig, f"{PDF_FIGS}/{img_path}.pdf", width=500, height=330)

    row, col = divmod(i, n_cols)

    for trace in subfig.data:
        fig.add_trace(trace, row=row + 1, col=col + 1)


# %%
# Update font size for each individual subplot title
for annotation in fig.layout.annotations:
    annotation.font.size = 12  # Adjust to your desired font size

fig.update_traces(showlegend=False)
for trace in fig.select_traces(row=1, col=1):
    trace.update(showlegend=True)
fig.update_layout(showlegend=True)

# place the legend above the subplots
fig.layout.legend.update(
    y=1.08, xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)", orientation="h"
)

# set the figure size based on the number of rows and columns
fig.layout.update(height=230 * n_rows)
fig.layout.update(width=280 * n_cols)

# set the shared y and x axis ranges to (-0.2, 0.2), and (0, 0.2)
fig.update_xaxes(range=[-0.2, 0.2])
fig.update_yaxes(range=[0, 0.2])

# Create shared x and y axis titles
if subfig is not None:
    x_title = subfig.layout.xaxis.title.text  # used in annotations below
    y_title = subfig.layout.yaxis.title.text
else:
    raise ValueError("x_title and y_title are not defined")

for i in range(1, n_rows + 1):
    for j in range(1, n_cols + 1):
        fig.update_xaxes(title_text="", row=i, col=j)
        fig.update_yaxes(title_text="", row=i, col=j)

axis_titles = dict(xref="paper", yref="paper", showarrow=False, font_size=16)
portrait = n_rows > n_cols
fig.add_annotation(  # x-axis title
    x=0.5,
    y=-0.09 if portrait else -0.18,
    text=x_title,
    borderpad=5,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=-0.09 if portrait else -0.07,
    y=0.5,
    text=y_title,
    textangle=-90,
    borderpad=5,
    **axis_titles,
)

# standardize the margins and template
fig.layout.margin.update(l=60, r=10, t=0 if portrait else 10, b=60 if portrait else 10)
fig.update_xaxes(matches=None)
fig.update_yaxes(matches=None)
fig.layout.template = "pymatviz_white"

fig.show()


# %%
img_path = (
    f"tile-rolling-mae-batches-{n_rows}x{n_cols}"
    f"{'-only-compliant' if not show_non_compliant else ''}"
)
pmv.save_fig(fig, f"{SITE_FIGS}/{img_path}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_path}.pdf")
