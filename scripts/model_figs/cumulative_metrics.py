"""Plot cumulative metrics like precision, recall, F1, MAE, RMSE as lines for all models
into face plot with one subplot per metric. Cumulative here means descending the list of
test set materials ranked by model-predicted stability starting from the most stable
and updating the metric (Recall, MAE, etc.) after each new material. This plot
simulates an actual materials screening process and allows practitioners to choose
a cutoff point for the number of DFT calculations they have budget and see which model
will provide the best hit rate for the given budget.
"""

# %%
import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import MATPLOTLIB, PLOTLY

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.plots import cumulative_metrics
from matbench_discovery.preds import df_each_pred, df_preds, model_styles, models

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %%
metrics: tuple[str, ...] = globals().get("metrics", ("Precision", "Recall"))
# metrics = ("MAE",)
range_y = {
    ("MAE",): (0, 0.7),
    ("Precision", "Recall"): (0, 1),
}[metrics]
fig, df_metric = cumulative_metrics(
    e_above_hull_true=df_preds[MbdKey.each_true],
    df_preds=df_each_pred[models],
    project_end_point="xy",
    backend=(backend := PLOTLY),
    metrics=metrics,
    # facet_col_wrap=2,
    # increase facet col gap
    facet_col_spacing=0.05,
    # markers=True,
    show_n_stable=metrics != ("MAE",),
)

x_label = "Number of screened materials"
if backend == MATPLOTLIB:
    # fig.suptitle(title)
    fig.text(0.5, -0.08, x_label, ha="center", fontdict={"size": 16})
if backend == PLOTLY:
    for key in filter(lambda key: key.startswith("yaxis"), fig.layout):
        fig.layout[key].range = range_y

    fig.layout.margin.update(l=0, r=0, t=30, b=50)
    # use annotation for x-axis label
    fig.add_annotation(
        **dict(x=0.5, y=-0.15, xref="paper", yref="paper"),
        text=x_label,
        showarrow=False,
        font=dict(size=14),
    )
    fig.update_traces(line=dict(width=3))
    fig.layout.legend.update(bgcolor="rgba(0,0,0,0)")
    # fig.layout.legend.update(
    #     orientation="h", yanchor="bottom", y=1.1, xanchor="center", x=0.5
    # )
    if "MAE" in metrics:
        fig.layout.legend.update(traceorder="reversed")
    if metrics == ("MAE",):
        fig.layout.legend.update(y=1, x=1, xanchor="right", yanchor="top")
    if len(metrics) * len(models) != len(fig.data):
        raise ValueError(
            "Expected one trace per model per metric, i.e. "
            f"{len(metrics) * len(models)} traces, got {len(fig.data)}"
        )

    for trace in fig.data:
        if line_style := model_styles.get(trace.name):
            ls, _marker, color = line_style
            trace.line = dict(color=color, dash=ls, width=2)

        # show only the N best models by default
        # if trace.name in df_metrics.T.sort_values("F1").index[:-6]:
        #     trace.visible = "legendonly"

        last_idx = pd.Series(trace.y).last_valid_index()
        last_x = trace.x[last_idx]
        last_y = trace.y[last_idx]
        color = dict(color=trace.line.color)
        subplot_col = int(trace.xaxis[-1].strip("x") or 1)

        # PRO can be rotated so text is parallel to line, more readable
        # CON remains visible when trace is hidden via legend
        # fig.add_annotation(
        #     x=last_x,
        #     y=last_y,
        #     text=trace.name,
        #     font=color,
        #     row=1,
        #     col=subplot_col,
        #     standoff=30,
        #     textangle=40 if subplot_col == 1 else -20,
        # )
        # PRO can be toggled with the line via legendgroup
        # CON unreadable due to lines overlapping
        # fig.add_scatter(
        #     x=[last_x],
        #     y=[last_y],
        #     mode="markers+text",
        #     text=trace.name,
        #     textposition="top center",
        #     textfont=color,
        #     marker=color,
        #     legendgroup=trace.name,
        #     showlegend=False,
        #     xaxis=trace.xaxis,  # add to the right subplot
        #     yaxis=trace.yaxis,
        #     hoverinfo="skip",
        # )

fig.show()


# %%
img_name = f"cumulative-{'-'.join(metrics).lower()}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=1000, height=400)
