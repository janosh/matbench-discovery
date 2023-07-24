"""Compile metrics and total run times for all models and export them to JSON, a
pandas-styled HTML table and a plotly figure.
"""


# %%
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier

from matbench_discovery import FIGS, PDF_FIGS
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.plots import df_to_pdf, df_to_svelte_table
from matbench_discovery.preds import df_metrics, df_metrics_10k, each_true_col

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %% add dummy classifier results to df_metrics
dummy_clf = DummyClassifier(strategy="stratified", random_state=0)
df_mp = pd.read_csv(DATA_FILES.mp_energies, index_col=0)
dummy_clf.fit(np.zeros_like(df_mp.energy_above_hull), df_mp.energy_above_hull == 0)
dummy_clf_preds = dummy_clf.predict(np.zeros(len(df_wbm)))
true_clf = df_wbm[each_true_col] < 0
each_true = df_wbm[each_true_col]

dummy_metrics = stable_metrics(
    each_true, np.array([1, -1])[dummy_clf_preds.astype(int)]
)

# important: regression metrics from dummy_clf are meaningless, we overwrite them with
# correct values here. don't remove!
dummy_metrics["DAF"] = 1
dummy_metrics["R2"] = 0
dummy_metrics["MAE"] = (each_true - each_true.mean()).abs().mean()
dummy_metrics["RMSE"] = ((each_true - each_true.mean()) ** 2).mean() ** 0.5

df_metrics["Dummy"] = dummy_metrics
# add 'global' to dummy to make explicit that this is not dummy on first 10k most stable
# predictions but on the whole dataset
df_metrics_10k["Dummy"] = dummy_metrics


# %%
ontology = {  # (training type, test type, model type)
    "ALIGNN": ("RS2RE", "IS2RE", "GNN"),
    "ALIGNN Pretrained": ("RS2RE", "IS2RE", "GNN"),
    "CHGNet": ("S2EFSM", "IS2RE-SR", "UIP-GNN"),
    "M3GNet": ("S2EFS", "IS2RE-SR", "UIP-GNN"),
    "MEGNet": ("RS2RE", "IS2E", "GNN"),
    "CGCNN": ("RS2RE", "IS2E", "GNN"),
    "CGCNN+P": ("S2RE", "IS2RE", "GNN"),
    "Wrenformer": ("RP2RE", "IP2E", "Transformer"),
    "BOWSR": ("RS2RE", "IS2RE-BO", "BO-GNN"),
    "Voronoi RF": ("RS2RE", "IS2E", "Fingerprint"),
    "M3GNet->MEGNet": ("S2EFS", "IS2RE-SR", "UIP-GNN"),
    "CHGNet->MEGNet": ("S2EFSM", "IS2RE-SR", "UIP-GNN"),
    "Dummy": ("", "", ""),
}
ontology_cols = ["Trained", "Deployed", "Model Class"]
df_ont = pd.DataFrame(ontology, index=ontology_cols)


# %%
R2_col = "R<sup>2</sup>"
higher_is_better = {*f"DAF {R2_col} Precision Recall F1 Accuracy TPR TNR TP TN".split()}
lower_is_better = {"MAE", "RMSE", "FPR", "FNR", "FP", "FN"}

# if True, make metrics-table-megnet-uip-combos.(svelte|pdf) for /si
make_uip_megnet_comparison = False
hide_metrics = "TP FN FP TN FNR FPR Recall Trained Deployed".split()

for label, df, extra_hide_metrics in (
    # hide redundant metrics (TPR = Recall, FPR = 1 - TNR, FNR = 1 - TPR)
    ("-first-10k", df_metrics_10k, ["TPR", "TNR"]),
    ("", df_metrics, []),
):
    df_table = pd.concat([df, df_ont]).rename(index={"R2": R2_col})
    df_table.index.name = "Model"

    drop_models = ["CHGNet + MEGNet", "M3GNet + MEGNet"]
    if make_uip_megnet_comparison:
        drop_models = [*{*df_table} - {*drop_models, "MEGNet", "M3GNet", "CHGNet"}]
        label += "-uip-megnet-combos"
        print(
            "hint: for make_uip_megnet_comparison, uncomment the lines chgnet_megnet "
            "and m3gnet_megnet in PredFiles"
        )
    df_filtered = df_table.T.drop(drop_models).drop(
        [*hide_metrics, *extra_hide_metrics],  # type: ignore
        axis="columns",
        errors="ignore",
    )
    styler = (
        df_filtered.style.format(
            # render integers without decimal places
            dict.fromkeys("TP FN FP TN".split(), "{:,.0f}"),
            precision=2,  # render floats with 2 decimals
            na_rep="",  # render NaNs as empty string
        ).background_gradient(
            cmap="viridis", subset=list(higher_is_better & {*df_filtered})
        )
        # reverse color map if lower=better
        .background_gradient(
            cmap="viridis_r", subset=list(lower_is_better & {*df_filtered})
        )
    )
    styles = {
        "": "font-family: sans-serif; border-collapse: collapse;",
        "td, th": "border: none; padding: 4px 6px; white-space: nowrap;",
        "th": "border: 1px solid; border-width: 1px 0; text-align: left;",
    }
    styler.set_table_styles([dict(selector=sel, props=styles[sel]) for sel in styles])
    styler.set_uuid("")

    #  export model metrics as styled HTML table and Svelte component
    # draw dotted line between classification and regression metrics
    df_to_svelte_table(
        styler,
        f"{FIGS}/metrics-table{label}.svelte",
        inline_props="class='roomy'",
        styles="#T_ :is(td, th):nth-last-child(4) { border-left: 1px dotted white; }",
    )
    try:
        df_to_pdf(styler, f"{PDF_FIGS}/metrics-table{label}.pdf")
    except ImportError as exc:
        print(f"df_to_pdf failed: {exc}")


# %%
# hide_rows = list(set(df_metrics) - set(df_metrics.T.F1.nlargest(6).index))
# styler.hide(hide_rows)  # show only the best models by F1 score
png_metrics = f"{PDF_FIGS}/metrics-table.png"
try:
    import dataframe_image

    dataframe_image.export(styler, png_metrics, dpi=300)
except ImportError:
    print("dataframe_image not installed, skipping png export")
