"""Check if points with large error compared to DFT but little disagreement between
models can pinpoint DFT calculation gone wrong.
"""

# %%
import pandas as pd
import plotly.express as px
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

test_subset = globals().get("test_subset", TestSubset.uniq_protos)


show_non_compliant = globals().get("show_non_compliant", cli_args.show_non_compliant)
models_to_plot = complete_models(show_non_compliant=show_non_compliant)
df_preds = load_df_wbm_with_preds(models=models_to_plot, subset=test_subset)
model_labels = [model.label for model in models_to_plot]


# %%
df_preds[MbdKey.each_mean_models] = (
    df_preds[model_labels].mean(axis=1)
    + df_preds[MbdKey.each_true]
    - df_preds[MbdKey.e_form_dft]
)
df_preds[MbdKey.model_std_each] = df_preds[model_labels].std(axis=1)

df_each_err = pd.DataFrame()
for model in model_labels:
    df_each_err[model] = df_preds[model] - df_preds[MbdKey.e_form_dft]

df_preds[MbdKey.each_err_models] = df_each_err.abs().mean(axis=1)
del df_each_err


# %% scatter plot of largest model errors vs. DFT hull distance
# while some points lie on a horizontal line of constant error, more follow the identity
# line showing models are biased to predict low energies likely as a result of training
# on MP which is highly low-energy enriched.
# also possible models failed to learn whatever physics makes these materials unstable
n_structs = 200
df_plot = df_preds.nlargest(n_structs, MbdKey.each_err_models).round(2)
y_max = df_plot[MbdKey.each_mean_models].max() + 0.6
y_min = df_plot[MbdKey.each_mean_models].min() - 0.3

fig = px.scatter(
    df_plot,
    x=MbdKey.each_true,
    y=MbdKey.each_mean_models,
    color=MbdKey.model_std_each,
    hover_name=Key.mat_id,
    hover_data=[Key.formula],
    color_continuous_scale="Turbo",
    range_y=(y_min, y_max),
)
# for horizontal colorbar
# yanchor="bottom", y=1, xanchor="center", x=0.5, orientation="h", thickness=12
fig.layout.coloraxis.colorbar.update(title_side="right", thickness=14)
fig.layout.margin.update(l=60, r=10, t=30, b=60)
pmv.powerups.add_identity_line(fig)

# size markers by square root of structure site count
fig.data[0].marker.size = df_plot["n_sites"] ** 0.5 * 3
fig.show()
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"scatter-largest-errors-models-mean-vs-true-hull-dist-all{img_suffix}"
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=600, height=300)
