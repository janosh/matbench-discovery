"""Check if points with large error compared to DFT but little disagreement between
models can pinpoint DFT calculation gone wrong.
"""


# %%
import pandas as pd
from crystal_toolkit.helpers.utils import hook_up_fig_with_struct_viewer
from pymatviz.utils import add_identity_line, save_fig

from matbench_discovery import FIGS, PDF_FIGS
from matbench_discovery.data import DATA_FILES
from matbench_discovery.preds import (
    df_preds,
    each_true_col,
    model_mean_each_col,
    model_mean_err_col,
    model_std_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"


# %% scatter plot of largest model errors vs. DFT hull distance
# while some points lie on a horizontal line of constant error, more follow the identity
# line showing models are biased to predict low energies likely as a result of training
# on MP which is highly low-energy enriched.
# also possible models failed to learn whatever physics makes these materials highly
# unstable

material_classes = {
    "all": r".*",
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
n_structs = 200

for material_cls, pattern in material_classes.items():
    df_subset = df_preds[df_preds["formula"].str.match(pattern)]
    df_plot = df_subset.nlargest(n_structs, model_mean_err_col).round(2)

    fig = df_plot.plot.scatter(
        x=each_true_col,
        y=model_mean_each_col,
        color=model_std_col,
        size="n_sites",
        backend="plotly",
        hover_name="material_id",
        hover_data=["formula"],
        color_continuous_scale="Turbo",
        range_x=[-0.5, 4],
        range_y=[-0.5, 4],
        # range_color=[0, df_plot[model_std_col].max()],
    )
    # for horizontal colorbar
    # yanchor="bottom", y=1, xanchor="center", x=0.5, orientation="h", thickness=12
    fig.layout.coloraxis.colorbar.update(title_side="right", thickness=14)
    fig.layout.margin.update(l=0, r=30, b=0, t=60)
    add_identity_line(fig)
    fig.layout.title.update(
        text=f"{n_structs} largest {material_cls} model errors: Predicted vs.<br>"
        "DFT hull distance colored by model disagreement",
        x=0.5,
    )
    # tried setting error_y=model_std_col but looks bad
    # fig.update_traces(
    #     error_y=dict(color="rgba(255,255,255,0.2)", width=3, thickness=2)
    # )
    fig.show()
    img_name = f"scatter-largest-errors-models-mean-vs-true-hull-dist-{material_cls}"
    save_fig(fig, f"{FIGS}/{img_name}.svelte")
    save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")


# %%
df_cse = pd.read_json(DATA_FILES.wbm_cses_plus_init_structs).set_index("material_id")


# %% struct viewer
app = hook_up_fig_with_struct_viewer(
    fig,
    df_cse,
    "initial_structure",
    # validate_id requires material_id to be hover_name
    validate_id=lambda id: id.startswith(("wbm-", "mp-", "mvc-")),
)
app.run(port=8000)
