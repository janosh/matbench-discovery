"""Plot cumulative precision and/or recall and/or F1 curves for all models into facet
plot with one subplot per metric. Cumulative here means going through the list of WBM
materials ranked by the model's stability prediction starting from the most stable
and updating the precision, recall and F1 score after each new material. This plot
simulates an actual materials screening process and allows practitioners to choose
a cutoff point for the number of DFT calculations they have budget and see which model
will provide the best hit rate for the given budget.
"""


# %%
import pandas as pd
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, ROOT
from matbench_discovery.plots import cumulative_precision_recall
from matbench_discovery.preds import df_each_pred, df_metrics, df_preds, each_true_col

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


# %%
fig, df_metric = cumulative_precision_recall(
    e_above_hull_true=df_preds[each_true_col],
    df_preds=df_each_pred,
    project_end_point="xy",
    backend=(backend := "plotly"),
    range_y=(0, 1)
    # template="plotly_white",
)

x_label = "Number of screened WBM materials"
if backend == "matplotlib":
    # fig.suptitle(title)
    fig.text(0.5, -0.08, x_label, ha="center", fontdict={"size": 16})
if backend == "plotly":
    fig.layout.legend.update(x=0, y=0, bgcolor="rgba(0,0,0,0)")
    fig.layout.margin.update(l=0, r=5, t=30, b=50)
    fig.add_annotation(
        x=0.5,
        y=-0.15,
        xref="paper",
        yref="paper",
        text=x_label,
        showarrow=False,
        font=dict(size=14),
    )
    fig.update_traces(line=dict(width=3))
    for trace in fig.data:
        if trace.name in df_metrics.T.sort_values("F1").index[:-6]:
            trace.visible = "legendonly"  # show only top models by default
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
img_name = "cumulative-clf-metrics"
save_fig(fig, f"{FIGS}/{img_name}.svelte")
# save_fig(fig, f"{STATIC}/{img_name}.webp", scale=3)
save_fig(fig, f"{ROOT}/paper/figures/{img_name}.pdf", width=720, height=370)
