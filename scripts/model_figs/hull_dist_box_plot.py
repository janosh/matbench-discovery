# %%
import plotly.express as px
import plotly.graph_objects as go
import pymatviz as pmv

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, Quantity, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.preds.discovery import df_each_err, df_preds, models

__author__ = "Janosh Riebesell"
__date__ = "2023-05-25"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
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

# Get the default Plotly colors that will be used for the boxes
color_seq = px.colors.qualitative.Plotly

for idx, model in enumerate(models_to_plot):
    ys = [df_each_err[model].quantile(quant) for quant in (0.05, 0.25, 0.5, 0.75, 0.95)]

    # Use the same color for both box and label
    color = color_seq[idx % len(color_seq)]
    fig.add_box(y=ys, name=model, width=0.8, marker_color=color)

    # annotate median with numeric value
    median = ys[2]
    fig.add_annotation(
        x=idx, y=median, text=f"{median:.2}", showarrow=False, font_size=9
    )

fig.layout.showlegend = False
# use line breaks to offset every other x-label and color them
x_labels_with_offset = [
    f"{'<br>' * (idx % 3)}<span style='color: {color_seq[idx % len(color_seq)]}'>"
    f"{label}</span>"
    for idx, label in enumerate(models_to_plot)
]

# prevent x-labels from rotating
fig.layout.xaxis.range = [-0.7, len(models_to_plot) - 0.3]
fig.layout.xaxis.update(
    tickangle=0,
    tickvals=models_to_plot,
    ticktext=x_labels_with_offset,
)
fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"box-hull-dist-errors{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
fig.layout.showlegend = False
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
