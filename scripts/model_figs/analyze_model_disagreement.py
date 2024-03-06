"""Check if points with large error compared to DFT but little disagreement between
models can pinpoint DFT calculation gone wrong.
"""

# %%
import sys

import pandas as pd
from pymatviz.io import save_fig
from pymatviz.utils import add_identity_line

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.data import DATA_FILES
from matbench_discovery.enums import Key, TestSubset
from matbench_discovery.preds import df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

test_subset = globals().get("test_subset", TestSubset.full)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)


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
n_structs, fig = 200, None

for material_cls, pattern in material_classes.items():
    df_subset = df_preds[df_preds[Key.formula].str.match(pattern)]
    df_plot = df_subset.nlargest(n_structs, Key.model_mean_err).round(2)

    fig = df_plot.plot.scatter(
        x=Key.each_true,
        y=Key.model_mean_each,
        color=Key.model_std_each,
        backend="plotly",
        hover_name=Key.mat_id,
        hover_data=[Key.formula],
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
    label = {"all": "structures"}.get(material_cls, material_cls)
    fig.layout.title.update(
        text=f"{n_structs} {material_cls} with largest hull distance errors<br>"
        "colored by model disagreement, sized by number of sites",
        x=0.5,
    )
    # size markers by structure
    fig.data[0].marker.size = df_plot["n_sites"] ** 0.5 * 3
    # tried setting error_y=model_std_col but looks bad
    # fig.update_traces(
    #     error_y=dict(color="rgba(255,255,255,0.2)", width=3, thickness=2)
    # )
    fig.show()
    img_name = f"scatter-largest-errors-models-mean-vs-true-hull-dist-{material_cls}"
    save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
    save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")


# %%
df_cse = pd.read_json(DATA_FILES.wbm_cses_plus_init_structs).set_index(Key.mat_id)


# %% CTK structure viewer
is_jupyter = "ipykernel" in sys.modules
if is_jupyter:  # only run this in Jupyter Notebook
    from crystal_toolkit.helpers.utils import hook_up_fig_with_struct_viewer

    app = hook_up_fig_with_struct_viewer(
        fig,
        df_cse,
        Key.init_struct,
        # validate_id requires material_id to be hover_name
        validate_id=lambda mat_id: mat_id.startswith(("wbm-", "mp-", "mvc-")),
    )
    app.run(port=8000)
