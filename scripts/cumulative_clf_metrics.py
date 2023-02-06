# %%
import pandas as pd
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, ROOT, STATIC
from matbench_discovery.plots import cumulative_precision_recall
from matbench_discovery.preds import df_each_pred, df_metrics, df_wbm, each_true_col

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


# %%
fig, df_metric = cumulative_precision_recall(
    e_above_hull_true=df_wbm[each_true_col],
    df_preds=df_each_pred[df_metrics.T.F1.nlargest(6).index],
    project_end_point="xy",
    backend=(backend := "plotly"),
    metrics=("Precision", "Recall"),
    # template="plotly_white",
)

# title = f"{today} - Cumulative Precision, Recall, F1 scores for classifying stable materials"
xlabel = "Number of WBM materials"
if backend == "matplotlib":
    # fig.suptitle(title)
    fig.text(0.5, -0.08, xlabel, ha="center", fontdict={"size": 16})
if backend == "plotly":
    fig.layout.legend.update(
        x=0.02, y=0.02, itemsizing="constant", bgcolor="rgba(0,0,0,0)"
    )  # , title=title
    # fig.layout.height = 500
    fig.layout.margin = dict(l=0, r=5, t=30, b=60)
    fig.add_annotation(
        x=0.5,
        y=-0.15,
        xref="paper",
        yref="paper",
        text=xlabel,
        showarrow=False,
        font=dict(size=16),
    )
    mode = "dark"
    fig.layout.template = f"plotly_{'white' if mode == 'light' else 'dark'}"
    fig.update_traces(line=dict(width=3))
    for trace in fig.data:
        last_idx = pd.Series(trace.y).last_valid_index()
        last_x = trace.x[last_idx]
        last_y = trace.y[last_idx]
        color = dict(color=trace.line.color)
        col = trace.xaxis[-1].strip("x") or 1

        # PRO can be rotated so text is parallel to line, more readable
        # CON remains visible when trace is hidden via legend
        # fig.add_annotation(
        #     x=last_x,
        #     y=last_y,
        #     text=trace.name,
        #     font=color,
        #     row=1,
        #     col=int(col),
        #     standoff=30,
        #     textangle=40 if col == 1 else -20,
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
# file will be served by site
# so we round y floats to reduce file size since
for trace in fig.data:
    assert isinstance(trace.y[0], float)
    trace.y = [round(y, 3) for y in trace.y]

img_name = "cumulative-clf-metrics"
save_fig(fig, f"{FIGS}/{img_name}.svelte")
save_fig(fig, f"{STATIC}/{img_name}.webp", scale=3)
save_fig(fig, f"{ROOT}/tmp/figures/{img_name}.pdf", width=700, height=350)
