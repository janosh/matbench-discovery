"""parity plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. First 2 plots put all models in single figure with selectable traces.
Last plot is split into 2x3 subplots, one for each model.
"""

# %%
import itertools
import os

import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2024-09-07"

# toggle between formation energy and energy above convex hull
update_existing: bool = cli_args.update_existing
models_to_plot = cli_args.models  # get model list from CLI, defaults to all models
energy_labels = {Key.e_form: "E<sub>form</sub>", Key.each: "E<sub>hull dist</sub>"}

# Load predictions for specified models
df_preds = load_df_wbm_with_preds(
    models=models_to_plot, subset=cli_args.test_subset
).round(3)


# %% parity plot of actual vs predicted e_form_per_atom
parity_scatters_dir = f"{SITE_FIGS}/energy-parity"
os.makedirs(parity_scatters_dir, exist_ok=True)

for model, which_energy in itertools.product(models_to_plot, (Key.e_form, Key.each)):
    energy_key = which_energy.replace("_", "-")
    img_name = f"{energy_key}-parity-{model.key.replace(' ', '-')}"
    img_path = f"{parity_scatters_dir}/{img_name}.svelte"
    if os.path.isfile(img_path) and not update_existing:
        print(f"{img_path} already exists, skipping")
        continue

    if which_energy == Key.each:
        # Calculate EACH prediction for this model
        each_pred = (
            df_preds[MbdKey.each_true]
            + df_preds[model.label]
            - df_preds[MbdKey.e_form_dft]
        )
        df_in = df_preds[[Key.formula, Key.mat_id, MbdKey.each_true]].copy()
        df_in[model.label] = each_pred
        e_true_col = MbdKey.each_true
    elif which_energy == Key.e_form:
        df_in = df_preds[[Key.formula, Key.mat_id, MbdKey.e_form_dft, model.label]]
        e_true_col = MbdKey.e_form_dft
    else:
        raise ValueError(f"Unexpected {which_energy=}")

    e_pred_col = f"{model.label} {e_true_col.label.replace('DFT ', '')}"
    df_in = df_in.rename(columns={model.label: e_pred_col})
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
    x_title = f"PBE {energy_labels[which_energy]} (eV/atom)"
    fig.layout.xaxis.title.update(text=x_title, font_size=16)
    y_title = f"{model.label} {energy_labels[which_energy]} (eV/atom)"
    fig.layout.yaxis.title.update(text=y_title, font_size=16)
    fig.layout.coloraxis.colorbar.update(
        title="Point Density", title_side="bottom", tickangle=0
    )

    fig.layout.title.update(text=f"{model.label} {which_energy}", x=0.5)
    fig.layout.margin.update(l=0, r=0, t=50, b=0)
    fig.show()
    pmv.save_fig(fig, img_path)
