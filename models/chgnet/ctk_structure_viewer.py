import pandas as pd
from crystal_toolkit.helpers.utils import hook_up_fig_with_struct_viewer
from pymatviz.enums import Key
from pymatviz.utils import PLOTLY

from matbench_discovery.data import Model

__author__ = "Janosh Riebesell"
__date__ = "2023-03-07"

"""
This scripts runs a Crystal Toolkit app that shows a parity plot of CHGNet vs PBE
energies and allows to click on points to view the corresponding structures. Run with:
python scripts/ctk_structure_viewer.py
Then open http://localhost:8000 in your browser.
"""

e_form_2000 = "e_form_per_atom_chgnet_2000"
e_form_500 = "e_form_per_atom_chgnet_500"

df_chgnet = pd.read_json(Model.chgnet.path.replace(".csv.gz", ".json.gz"))
df_chgnet = df_chgnet.set_index(Key.mat_id)

df_chgnet_2000 = pd.read_csv(Model.chgnet.path)
df_chgnet_2000 = df_chgnet_2000.set_index(Key.mat_id).add_suffix("_2000")
df_chgnet[list(df_chgnet_2000)] = df_chgnet_2000

df_chgnet_500 = pd.read_csv(Model.chgnet.path.replace("-06", "-04"))
df_chgnet_500 = df_chgnet_500.set_index(Key.mat_id).add_suffix("_500")
df_chgnet[list(df_chgnet_500)] = df_chgnet_500

e_form_abs_diff = "e_form_abs_diff"
min_e_diff = 0.1
df_chgnet[e_form_abs_diff] = abs(df_chgnet[e_form_2000] - df_chgnet[e_form_500])
df_plot = df_chgnet.round(3).query(f"{e_form_abs_diff} > {min_e_diff}")


plot_labels = {
    e_form_500: "CHGNet E<sub>form</sub> after 500 steps",
    e_form_2000: "CHGNet E<sub>form</sub> after 2000 steps",
    e_form_abs_diff: "Δ E<sub>form</sub>",
}

fig = df_plot.reset_index().plot.scatter(
    x=e_form_500,
    y=e_form_2000,
    backend=PLOTLY,
    hover_name=Key.mat_id,
    hover_data=[Key.formula],
    labels=plot_labels,
    size=e_form_abs_diff,
    color=e_form_abs_diff,
    template="plotly_white",
)

fig.layout.margin.update(b=20, l=40, r=20, t=50)
fig.layout.coloraxis.colorbar.update(
    title=dict(text="Energy Diff (eV/atom)", side="right"), thickness=10
)
# slightly increase scatter point size (lower sizeref means larger)
fig.update_traces(marker_sizeref=0.02, selector=dict(mode="markers"))

app = hook_up_fig_with_struct_viewer(
    fig,
    df_plot,
    "chgnet_structure",
    # validate_id requires material_id to be hover_name
    validate_id=lambda mat_id: mat_id.startswith(("wbm-", "mp-", "mvc-")),
)
app.run()
