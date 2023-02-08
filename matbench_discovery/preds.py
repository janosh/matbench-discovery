from __future__ import annotations

import pandas as pd

from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_preds
from matbench_discovery.metrics import stable_metrics

"""Centralize data-loading and computing metrics for plotting scripts"""

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"

models = sorted(
    "Wrenformer, CGCNN+P, Voronoi Random Forest, MEGNet, M3GNet + MEGNet, "
    "BOWSR + MEGNet".split(", ")
)
e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
each_pred_col = "e_above_hull_pred"

df_wbm = load_df_wbm_preds(list(PRED_FILENAMES)).round(3)
drop_cols = {*PRED_FILENAMES} - {*models}


df_metrics = pd.DataFrame()
df_metrics.index.name = "model"
for model in list(PRED_FILENAMES):
    df_metrics[model] = stable_metrics(
        df_wbm[each_true_col],
        df_wbm[each_true_col] + df_wbm[model] - df_wbm[e_form_col],
    )

# pick F1 as primary metric to sort by
df_metrics = df_metrics.round(3).sort_values("F1", axis=1)


df_each_pred = pd.DataFrame()
for model in df_metrics.T.MAE.sort_values().index:
    df_each_pred[model] = df_wbm[each_true_col] + df_wbm[model] - df_wbm[e_form_col]
