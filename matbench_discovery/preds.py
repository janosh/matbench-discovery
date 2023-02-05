from __future__ import annotations

import pandas as pd

from matbench_discovery.data import load_df_wbm_preds
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

df_wbm = load_df_wbm_preds(models).round(3)

for col in [e_form_col, each_true_col]:
    assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


df_metrics = pd.DataFrame()
for model in models:
    df_metrics[model] = stable_metrics(
        df_wbm[each_true_col],
        df_wbm[each_true_col] + df_wbm[model] - df_wbm[e_form_col],
    )

assert df_metrics.T.MAE.between(0, 0.2).all(), "MAE not in range"
assert df_metrics.T.R2.between(0.1, 1).all(), "R2 not in range"
assert df_metrics.T.RMSE.between(0, 0.25).all(), "RMSE not in range"
assert df_metrics.isna().sum().sum() == 0, "NaNs in metrics"
