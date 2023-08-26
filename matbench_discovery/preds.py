from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pandas as pd
from tqdm import tqdm

from matbench_discovery import ROOT, STABILITY_THRESHOLD
from matbench_discovery.data import Files, df_wbm, glob_to_df
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.plots import (
    ev_per_atom,
    model_labels,
    plotly_colors,
    plotly_line_styles,
    plotly_markers,
    quantity_labels,
)

"""Centralize data-loading and computing metrics for plotting scripts"""

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"

e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
each_pred_col = "e_above_hull_pred"
model_mean_each_col = "Mean prediction all models"
model_mean_err_col = "Mean error all models"
model_std_col = "Std. dev. over models"

for col in (model_mean_each_col, model_mean_err_col, model_std_col):
    quantity_labels[col] = f"{col} {ev_per_atom}"


class PredFiles(Files):
    """Data files provided by Matbench Discovery.
    See https://janosh.github.io/matbench-discovery/contribute for data descriptions.
    """

    # BOWSR optimizer coupled with original megnet
    bowsr_megnet = "bowsr/2023-01-23-bowsr-megnet-wbm-IS2RE.csv.gz"
    # default CHGNet model from publication with 400,438 params
    chgnet = "chgnet/2023-03-06-chgnet-wbm-IS2RE.csv.gz"

    # CGCnn 10-member ensemble
    cgcnn = "cgcnn/2023-01-26-cgcnn-ens=10-wbm-IS2RE.csv.gz"
    # CGCnn 10-member ensemble with 5-fold training set perturbations
    cgcnn_p = "cgcnn/2023-02-05-cgcnn-perturb=5-wbm-IS2RE.csv.gz"

    # original M3GNet straight from publication, not re-trained
    m3gnet = "m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv.gz"
    # m3gnet_direct = "m3gnet/2023-05-30-m3gnet-direct-wbm-IS2RE.csv.gz"
    # m3gnet_ms = "m3gnet/2023-06-01-m3gnet-manual-sampling-wbm-IS2RE.csv.gz"

    # MACE trained on original M3GNet training set
    mace = "mace/2023-08-14-mace-wbm-IS2RE-FIRE.csv.gz"

    # original MEGNet straight from publication, not re-trained
    megnet = "megnet/2022-11-18-megnet-wbm-IS2RE.csv.gz"
    # # CHGNet-relaxed structures fed into MEGNet for formation energy prediction
    # chgnet_megnet = "chgnet/2023-03-04-chgnet-wbm-IS2RE.csv.gz"
    # # M3GNet-relaxed structures fed into MEGNet for formation energy prediction
    # m3gnet_megnet = "m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv.gz"
    # megnet_rs2re = "megnet/2023-08-23-megnet-wbm-RS2RE.csv.gz"

    # Magpie composition+Voronoi tessellation structure features + sklearn random forest
    voronoi_rf = "voronoi/2022-11-27-train-test/e-form-preds-IS2RE.csv.gz"

    # wrenformer 10-member ensemble
    wrenformer = "wrenformer/2022-11-15-wrenformer-ens=10-IS2RE-preds.csv.gz"

    alignn = "alignn/2023-06-02-alignn-wbm-IS2RE.csv.gz"
    # alignn_pretrained = "alignn/2023-06-03-mp-e-form-alignn-wbm-IS2RE.csv.gz"
    # alignn_ff = "alignn_ff/2023-07-11-alignn-ff-wbm-IS2RE.csv.gz"


# model_labels remaps model keys to pretty plot labels (see Files)
PRED_FILES = PredFiles(root=f"{ROOT}/models", key_map=model_labels)


def load_df_wbm_with_preds(
    models: Sequence[str] = (*PRED_FILES,),
    pbar: bool = True,
    id_col: str = "material_id",
    **kwargs: Any,
) -> pd.DataFrame:
    """Load WBM summary dataframe with model predictions from disk.

    Args:
        models (Sequence[str], optional): Model names must be keys of
            matbench_discovery.data.PRED_FILES. Defaults to all models.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        id_col (str, optional): Column to set as df.index. Defaults to "material_id".
        **kwargs: Keyword arguments passed to glob_to_df().

    Raises:
        ValueError: On unknown model names.

    Returns:
        pd.DataFrame: WBM summary dataframe with model predictions.
    """
    if mismatch := ", ".join(set(models) - set(PRED_FILES)):
        raise ValueError(
            f"Unknown models: {mismatch}, expected subset of {set(PRED_FILES)}"
        )

    dfs: dict[str, pd.DataFrame] = {}

    try:
        for model_name in (bar := tqdm(models, disable=not pbar, desc="Loading preds")):
            bar.set_postfix_str(model_name)
            df = glob_to_df(PRED_FILES[model_name], pbar=False, **kwargs)
            df = df.set_index(id_col)
            dfs[model_name] = df
    except Exception as exc:
        raise RuntimeError(f"Failed to load {model_name=}") from exc

    from matbench_discovery.data import df_wbm

    df_out = df_wbm.copy()
    for model_name, df in dfs.items():
        model_key = model_name.lower().replace("→", "_").replace(" ", "_")

        cols = [col for col in df if col.startswith(f"e_form_per_atom_{model_key}")]
        if cols:
            if len(cols) > 1:
                print(
                    f"Warning: multiple pred cols for {model_name=}, using {cols[0]!r} "
                    f"out of {cols=}"
                )
            df_out[model_name] = df[cols[0]]

        elif pred_cols := list(df.filter(like="_pred_ens")):
            assert len(pred_cols) == 1
            df_out[model_name] = df[pred_cols[0]]
            if std_cols := list(df.filter(like="_std_ens")):
                df_out[f"{model_name}_std"] = df[std_cols[0]]

        elif pred_cols := list(df.filter(like=r"_pred_")):
            # make sure we average the expected number of ensemble member predictions
            assert len(pred_cols) == 10, f"{len(pred_cols) = }, expected 10"
            df_out[model_name] = df[pred_cols].mean(axis=1)

        else:
            cols = list(df)
            raise ValueError(
                f"No pred col for {model_name=} ({model_key=}), available {cols=}"
            )

    return df_out


# load WBM summary dataframe with all models' formation energy predictions (eV/atom)
df_preds = load_df_wbm_with_preds().round(3)
# for combo in [["CHGNet", "M3GNet"]]:
#     df_preds[" + ".join(combo)] = df_preds[combo].mean(axis=1)
#     PRED_FILES[" + ".join(combo)] = "combo"


df_metrics = pd.DataFrame()
df_metrics_10k = pd.DataFrame()  # look only at each model's 10k most stable predictions
prevalence = (df_wbm[each_true_col] <= STABILITY_THRESHOLD).mean()

df_metrics.index.name = "model"
for model in PRED_FILES:
    each_pred = df_preds[each_true_col] + df_preds[model] - df_preds[e_form_col]
    df_metrics[model] = stable_metrics(df_preds[each_true_col], each_pred)
    most_stable_10k = each_pred.nsmallest(10_000)
    df_metrics_10k[model] = stable_metrics(
        df_preds[each_true_col].loc[most_stable_10k.index], most_stable_10k
    )
    df_metrics_10k[model]["DAF"] = df_metrics_10k[model]["Precision"] / prevalence


# pick F1 as primary metric to sort by
df_metrics = df_metrics.round(3).sort_values("F1", axis=1, ascending=False)
df_metrics_10k = df_metrics_10k.round(3).sort_values("F1", axis=1, ascending=False)

models = list(df_metrics.T.MAE.sort_values().index)
# used for consistent markers, line styles and colors for a given model across plots
model_styles = dict(zip(models, zip(plotly_line_styles, plotly_markers, plotly_colors)))

# To avoid confusion for anyone reading this code, we calculate the formation energy MAE
# here and report it as the MAE for the energy above the convex hull prediction. The
# former is more easily calculated but the two quantities are the same. The formation
# energy of a material is the difference in energy between a material and its
# constituent elements in their standard states. The distance to the convex hull is
# defined as the difference between a material's formation energy and the minimum
# formation energy of all possible stable materials made from the same elements. Since
# the formation energy of a material is used to calculate the distance to the convex
# hull, the error of a formation energy prediction directly determines the error in the
# distance to the convex hull prediction.

# A further point of clarification: whenever we say convex hull distance we mean
# the signed distance that is positive for thermodynamically unstable materials above
# the hull and negative for stable materials below it.

# dataframe of all models' energy above convex hull (EACH) predictions (eV/atom)
df_each_pred = pd.DataFrame()
for model in models:
    df_each_pred[model] = (
        df_preds[each_true_col] + df_preds[model] - df_preds[e_form_col]
    )

# important: do df_each_pred.std(axis=1) before inserting
# df_each_pred[model_mean_each_col]
df_preds[model_std_col] = df_each_pred.std(axis=1)
df_each_pred[model_mean_each_col] = df_preds[model_mean_each_col] = df_each_pred.mean(
    axis=1
)

# dataframe of all models' errors in their EACH predictions (eV/atom)
df_each_err = pd.DataFrame()
for model in models:
    df_each_err[model] = df_preds[model] - df_preds[e_form_col]

df_each_err[model_mean_err_col] = df_preds[model_mean_err_col] = df_each_err.abs().mean(
    axis=1
)
