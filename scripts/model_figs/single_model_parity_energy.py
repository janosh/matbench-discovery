"""parity plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. First 2 plots put all models in single figure with selectable traces.
Last plot is split into 2x3 subplots, one for each model.
"""

# %%
import itertools
import os
from typing import Literal, get_args

import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import SITE_FIGS
from matbench_discovery.data import Model
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import (
    df_each_pred,
    df_metrics,
    df_metrics_uniq_protos,
    df_preds,
)

__author__ = "Janosh Riebesell"
__date__ = "2024-09-07"

# toggle between formation energy and energy above convex hull
EnergyType = Literal["e-form", "each"]
use_e_form, use_each = get_args(EnergyType)
energy_labels = {use_e_form: "Formation Energy", use_each: "Convex Hull Distance"}
update_existing: bool = False

test_subset = globals().get("test_subset", TestSubset.uniq_protos)
if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_metrics = df_metrics_uniq_protos


# %% parity plot of actual vs predicted e_form_per_atom
parity_scatters_dir = f"{SITE_FIGS}/energy-parity"
os.makedirs(parity_scatters_dir, exist_ok=True)

for model_name, which_energy in itertools.product(df_metrics, (use_e_form, use_each)):
    model_key = Model.from_label(model_name).key
    img_name = f"{which_energy}-parity-{model_name.lower().replace(' ', '-')}"
    img_path = f"{parity_scatters_dir}/{img_name}.svelte"
    if os.path.isfile(img_path) and not update_existing:
        continue

    if which_energy == use_each:
        df_in = df_each_pred.copy()
        df_in[MbdKey.each_true] = df_preds[MbdKey.each_true]
        df_in[Key.formula] = df_preds[Key.formula]
        df_in[Key.mat_id] = df_preds[Key.mat_id]
        e_true_col = MbdKey.each_true
    elif which_energy == use_e_form:
        df_in = df_preds
        e_true_col = MbdKey.e_form_dft
    else:
        raise ValueError(f"Unexpected {which_energy=}")

    e_pred_col = f"{model_name} {e_true_col.label.replace('DFT ', '')}"
    df_in = df_in.rename(columns={model_name: e_pred_col})
    n_points = len(df_in.dropna(subset=[e_true_col, e_pred_col]))

    fig = pmv.density_scatter_plotly(
        df=df_in.reset_index(drop=True),
        x=e_true_col,
        y=e_pred_col,
        hover_data=[Key.formula],
        hover_name=Key.mat_id,
        opacity=0.7,
        color_continuous_scale="agsunset",
        colorbar_kwargs=dict(orientation="h", thickness=15, x=0.3, y=0.8, len=0.5),
        stats=dict(
            prefix=f"N={n_points:,}<br>",
            x=0.99,
            xanchor="right",
            y=0.07,
            font_color="black",
        ),
        best_fit_line=dict(annotate_params=dict(y=0.01, font_size=16, x=0.99)),
    )
    fig.layout.xaxis.title.update(
        text=f"PBE {energy_labels[which_energy]} (eV/atom)", font_size=16
    )
    fig.layout.yaxis.title.update(
        text=f"{model_name} {energy_labels[which_energy]} (eV/atom)", font_size=16
    )

    fig.layout.coloraxis.colorbar.update(
        title="Point Density", title_side="bottom", tickangle=0
    )
    pmv.powerups.add_identity_line(fig)
    # fig.layout.update(width=600, height=400)

    pmv.save_fig(fig, img_path)
    fig.layout.title.update(text=f"{model_name} {which_energy}", x=0.5)
    fig.layout.margin.update(l=0, r=0, t=50, b=0)
    fig.show()
