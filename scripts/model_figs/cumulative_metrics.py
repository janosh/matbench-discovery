"""Plot cumulative metrics like precision, recall, F1, MAE, RMSE as lines for all models
into face plot with one subplot per metric. Cumulative here means descending the list of
test set materials ranked by model-predicted stability starting from the most stable
and updating the metric (Recall, MAE, etc.) after each new material. This plot
simulates an actual materials screening process and allows practitioners to choose
a cutoff point for the number of DFT calculations they have budget and see which model
will provide the best hit rate for the given budget.
"""

# %%
import numpy as np
import pymatviz as pmv
import scipy.interpolate

from matbench_discovery import PDF_FIGS, STABILITY_THRESHOLD, figs
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics.discovery import classify_stable
from matbench_discovery.plots import cumulative_metrics
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %%
metrics: tuple[str] | tuple[str, str] = globals().get(
    "metrics", ("Precision", "Recall")
)
# metrics = ("MAE",)
range_y = {
    ("MAE",): (0, 0.7),
    ("Precision", "Recall"): (0, 1),
}[metrics]

show_non_compliant = globals().get("show_non_compliant", cli_args.show_non_compliant)
models_to_plot = [
    model.label for model in complete_models(show_non_compliant=show_non_compliant)
]

fig, _df_metric = cumulative_metrics(
    e_above_hull_true=df_preds[MbdKey.each_true],
    df_preds=df_each_pred[models_to_plot],
    metrics=metrics,
    facet_col_spacing=0.05,  # increase facet col gap
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
if metrics == ("Precision", "Recall") and show_non_compliant:
    # site payload = full model set. curves are recomputed on a model-intrinsic grid
    # (not extracted from the figure, whose shared x grid depends on which models are
    # in the run) so an entry is identical in full regens and single-model merge runs
    cum_pr_models = []
    for label in models_to_plot:
        each_pred = df_each_pred[label].sort_values()
        each_true = df_preds[MbdKey.each_true].loc[each_pred.index]
        true_pos, false_neg, false_pos, _true_neg = classify_stable(
            each_true, each_pred, stability_threshold=STABILITY_THRESHOLD
        )
        n_true_pos_cum = true_pos.cumsum()  # all pd.Series, cumsum stays a Series
        precision_cum = n_true_pos_cum / (n_true_pos_cum + false_pos.cumsum())
        recall_cum = n_true_pos_cum / (n_true_pos_cum + false_neg.cumsum()).iloc[-1]
        # number of materials the model predicts stable = where its curve ends
        n_pred_stable = int((each_pred <= STABILITY_THRESHOLD).sum())
        if n_pred_stable < 2:  # can't happen for real models (thousands stable)
            raise ValueError(f"{label} predicts {n_pred_stable} stable materials")
        # log2-spaced sampling for higher density at the start of the discovery
        # campaign where metrics fluctuate most (mirrors plots.cumulative_metrics).
        # rounded to ints since x counts screened materials (also compresses better)
        log_xs = np.logspace(0, np.log2(n_pred_stable - 1), 100, base=2)
        xs = np.unique([*log_xs.round().astype(int), n_pred_stable])
        model_range = np.arange(n_pred_stable) + 1
        spline_degree = min(3, n_pred_stable - 1)  # k must be < n curve points
        precision, recall = (
            scipy.interpolate.make_interp_spline(
                model_range, curve.to_numpy()[:n_pred_stable], k=spline_degree
            )(xs)
            for curve in (precision_cum, recall_cum)
        )
        cum_pr_models.append(
            {
                "key": Model.from_label(label).key,
                "label": label,
                "x": figs.round_list(xs),
                "precision": figs.round_list(precision),
                "recall": figs.round_list(recall),
                # [n materials predicted stable, precision there, recall there]
                "end": [
                    n_pred_stable,
                    round(float(precision_cum.iloc[n_pred_stable - 1]), 5),
                    round(float(recall_cum.iloc[n_pred_stable - 1]), 5),
                ],
            }
        )
    n_stable = int((df_preds[MbdKey.each_true] <= STABILITY_THRESHOLD).sum())
    figs.write_site_payload(
        "cumulative-precision-recall",
        {"n_stable": n_stable, "models": cum_pr_models},
        assign_colors=True,
    )
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
