"""Plot cumulative metrics like precision, recall, F1, MAE, RMSE as lines for all models
into face plot with one subplot per metric. Cumulative here means descending the list of
test set materials ranked by model-predicted stability starting from the most stable
and updating the metric (Recall, MAE, etc.) after each new material. This plot
simulates an actual materials screening process and allows practitioners to choose
a cutoff point for the number of DFT calculations they have budget and see which model
will provide the best hit rate for the given budget.
"""

# %%
import pymatviz as pmv

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.plots import cumulative_metrics
from matbench_discovery.preds.discovery import df_each_pred, df_preds, models

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %%
metrics: tuple[str, ...] = globals().get("metrics", ("Precision", "Recall"))
# metrics = ("MAE",)
range_y = {
    ("MAE",): (0, 0.7),
    ("Precision", "Recall"): (0, 1),
}[metrics]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]

fig, df_metric = cumulative_metrics(
    e_above_hull_true=df_preds[MbdKey.each_true],
    df_preds=df_each_pred[models_to_plot],
    metrics=metrics,
    # facet_col_wrap=2,
    # increase facet col gap
    facet_col_spacing=0.05,
    # markers=True,
    endpoint_markers=(endpoint_markers := True),
    show_n_stable=metrics != ("MAE",),
)

x_label = "Number of screened materials"

for key in filter(lambda key: key.startswith("yaxis"), fig.layout):
    fig.layout[key].range = range_y

fig.layout.margin.update(l=0, r=0, t=20, b=0)
# use annotation for x-axis label
fig.add_annotation(
    **dict(x=0.5, y=-0.15, xref="paper", yref="paper"),
    text=x_label,
    showarrow=False,
    font=dict(size=14),
)
fig.update_traces(line=dict(width=3))
fig.layout.legend.update(bgcolor="rgba(0,0,0,0)")

if metrics == ("MAE",):
    fig.layout.legend.update(traceorder="reversed")
    # fig.layout.legend.update(y=1, x=1, xanchor="right", yanchor="top")

if len(metrics) * len(models_to_plot) * (2 if endpoint_markers else 1) != len(fig.data):
    raise ValueError(
        "Expected one trace per model per metric, i.e. "
        f"{len(metrics) * len(models_to_plot)} traces, got {len(fig.data)}"
    )

fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"cumulative-{'-'.join(metrics).lower()}{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
