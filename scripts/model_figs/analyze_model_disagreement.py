"""Check if points with large error compared to DFT but little disagreement between
models can pinpoint DFT calculation gone wrong.
"""

# %%
import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

test_subset = globals().get("test_subset", TestSubset.uniq_protos)


show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in cli_args.models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model.label])
]
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
material_classes = {"all": r".*"}
if split_out_classes := False:
    material_classes |= {
        "lanthanides": r".*(La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu)\d.*",
        "actinides": r".*(Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr)\d.*",
        "oxides": r".*O\d.*",
        "nitrides": r".*N\d.*",
        "sulfides": r".*S\d.*",
        "halides": r".*[FClBrI]\d.*",
        "pnictides": r".*[AsSbBi]\d.*",
        "chalcogenides": r".*[SeTe]\d.*",
        "borides": r".*B\d.*",
        "carbides": r".*C\d.*",
        "hydrides": r".*H\d.*",
        "oxynitrides": r".*[ON]\d.*",
    }
n_structs, fig = 200, None

for material_cls, pattern in material_classes.items():
    df_subset = df_preds[df_preds[Key.formula].str.match(pattern)]
    df_plot = df_subset.nlargest(n_structs, MbdKey.each_err_models).round(2)
    y_max = df_plot[MbdKey.each_mean_models].max() + 0.6
    y_min = df_plot[MbdKey.each_mean_models].min() - 0.3

    fig = df_plot.plot.scatter(
        x=MbdKey.each_true,
        y=MbdKey.each_mean_models,
        color=MbdKey.model_std_each,
        backend="plotly",
        hover_name=Key.mat_id,
        hover_data=[Key.formula],
        color_continuous_scale="Turbo",
        range_y=(-0.5, y_max),
    )
    # for horizontal colorbar
    # yanchor="bottom", y=1, xanchor="center", x=0.5, orientation="h", thickness=12
    fig.layout.coloraxis.colorbar.update(title_side="right", thickness=14)
    fig.layout.margin.update(l=60, r=10, t=30, b=60)
    pmv.powerups.add_identity_line(fig)
    label = {"all": "structures"}.get(material_cls, material_cls)
    # title = (
    #     f"{n_structs} {material_cls} with largest hull distance errors<br>"
    #     "colored by model disagreement, sized by number of sites"
    # )
    # fig.layout.title.update(text=title, x=0.5)

    # size markers by square root of structure site count
    fig.data[0].marker.size = df_plot["n_sites"] ** 0.5 * 3
    # tried setting error_y=model_std_col but looks bad
    # fig.update_traces(
    #     error_y=dict(color="rgba(255,255,255,0.2)", width=3, thickness=2)
    # )
    fig.show()
    img_suffix = "" if show_non_compliant else "-only-compliant"
    img_name = f"scatter-largest-errors-models-mean-vs-true-hull-dist-{material_cls}{img_suffix}"  # noqa: E501
    pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
    pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=600, height=300)
