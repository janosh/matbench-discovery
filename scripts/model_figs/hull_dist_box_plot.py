# %%
import plotly.graph_objects as go
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import Quantity, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.preds.discovery import df_each_err, df_preds, models

__author__ = "Janosh Riebesell"
__date__ = "2023-05-25"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_err = df_each_err.loc[df_preds.index]


# %%
show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]

fig = go.Figure()
fig.layout.yaxis.title = Quantity.e_above_hull_error
fig.layout.margin = dict(l=0, r=0, b=0, t=0)

for idx, model in enumerate(models_to_plot):
    ys = [df_each_err[model].quantile(quant) for quant in (0.05, 0.25, 0.5, 0.75, 0.95)]

    fig.add_box(y=ys, name=model, width=0.8)

    # annotate median with numeric value
    median = ys[2]
    fig.add_annotation(
        x=idx,
        y=median,
        text=f"{median:.2}",
        showarrow=False,
        # bgcolor="rgba(0, 0, 0, 0.2)",
    )

fig.layout.showlegend = False
# use line breaks to offset every other x-label
x_labels_with_offset = [
    f"{'<br>' * (idx % 2)}{label}" for idx, label in enumerate(models_to_plot)
]
# prevent x-labels from rotating
fig.layout.xaxis.range = [-0.7, len(models_to_plot) - 0.3]
fig.layout.xaxis.update(
    tickangle=0, tickvals=models_to_plot, ticktext=x_labels_with_offset
)
fig.layout.width = 70 * len(models_to_plot)
fig.show()


# %%
img_name = f"box-hull-dist-errors{'-only-compliant' if not show_non_compliant else ''}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
fig.layout.showlegend = False
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
