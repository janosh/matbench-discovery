from __future__ import annotations

import pandas as pd

from matbench_discovery.data import PRED_FILES, load_df_wbm_preds
from matbench_discovery.metrics import stable_metrics

"""Centralize data-loading and computing metrics for plotting scripts"""

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"

e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
each_pred_col = "e_above_hull_pred"

df_wbm = load_df_wbm_preds().round(3)


df_metrics = pd.DataFrame()
df_metrics.index.name = "model"
for model in list(PRED_FILES):
    df_metrics[model] = stable_metrics(
        df_wbm[each_true_col],
        df_wbm[each_true_col] + df_wbm[model] - df_wbm[e_form_col],
    )

# pick F1 as primary metric to sort by
df_metrics = df_metrics.round(3).sort_values("F1", axis=1)

# dataframe of all models' energy above convex hull (EACH) predictions (eV/atom)
df_each_pred = pd.DataFrame()
for model in df_metrics.T.MAE.sort_values().index:
    df_each_pred[model] = df_wbm[each_true_col] + df_wbm[model] - df_wbm[e_form_col]


# dataframe of all models' errors in their EACH predictions (eV/atom)
df_each_err = pd.DataFrame()
for model in df_metrics.T.MAE.sort_values().index:
    df_each_err[model] = df_wbm[model] - df_wbm[e_form_col]
