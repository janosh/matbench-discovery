# %%
import os
from collections.abc import Callable

import crystal_toolkit.components as ctc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from chgnet.model import StructOptimizer as ChgnetRelaxer
from chgnet.model.dynamics import TrajectoryObserver
from crystal_toolkit.settings import SETTINGS
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from m3gnet.models import Relaxer as M3gnetRelaxer
from pymatgen.core import Lattice, Structure
from pymatviz import IS_IPYTHON
from pymatviz.enums import Key

from matbench_discovery.data import df_wbm

# ruff: noqa: T201
__author__ = "Janosh Riebesell"
__date__ = "2023-03-22"

module_dir = os.path.dirname(__file__)


# %%
df_cse = pd.read_json(f"{module_dir}/wbm-chgnet-bad-relax.json.gz").set_index(
    Key.mat_id
)


# %%
material_id = "wbm-1-7168"
init_struct = Structure.from_dict(df_cse[Key.init_struct].loc[material_id])

final_struct = Structure.from_dict(df_cse[Key.cse].loc[material_id]["structure"])

init_spg = init_struct.get_space_group_info()
final_spg = final_struct.get_space_group_info()
print(f"{init_spg=}\n{final_spg=}")


# %%
chgnet_traj = ChgnetRelaxer().relax(init_struct)["trajectory"]
m3gnet_traj = M3gnetRelaxer().relax(init_struct)["trajectory"]


# %%
e_col = "Energy (eV)"
force_col = "Force (meV/Å)"
vol_col = "Volume (Å<sup>3</sup)"
df_chgnet, df_m3gnet = pd.DataFrame(), pd.DataFrame()

for df, traj in ((df_chgnet, chgnet_traj), (df_m3gnet, m3gnet_traj)):
    # convert array[float] to single float (only needed for m3gnet)
    df[e_col] = list(map(float, traj.energies))
    df[force_col] = [
        1000
        * np.linalg.norm(force, axis=1).mean()  # mean of norm of force on each atom
        for force in traj.forces
    ]
    df[vol_col] = [Lattice(cell).volume for cell in traj.cells]
    df.index.name = "Step"


# %%
row = df_wbm.loc[material_id]
dft_energy = row.uncorrected_energy + row.e_correction_per_atom_mp2020 * row.n_sites
print(f"{dft_energy=}")

# don't rerun this cell (see duplicate id error below)
struct_layouts: dict[str, html.Div] = {}
struct_comps: dict[str, ctc.StructureMoleculeComponent] = {}


# %%
def plot_energy_and_forces(
    df: pd.DataFrame,
    step: int,
    e_col: str,
    force_col: str,
    title: str,
) -> go.Figure:
    """Plot energy and forces as a function of relaxation step."""
    fig = go.Figure()
    # energy trace = primary y-axis
    fig.add_scatter(x=df.index, y=df[e_col], mode="lines", name="Energy")
    # get energy line color
    line_color = fig.data[0].line.color

    # forces trace = secondary y-axis
    fig.add_scatter(
        x=df.index,
        y=df[force_col],
        mode="lines",
        name="Forces",
        yaxis="y2",
    )

    fig.update_layout(
        template="plotly_white",
        title=title,
        xaxis={"title": "Relaxation Step"},
        yaxis={"title": e_col},
        yaxis2={"title": force_col, "overlaying": "y", "side": "right"},
        legend=dict(yanchor="top", y=1, xanchor="right", x=1),
    )

    # vertical line at the specified step
    fig.add_vline(x=step, line={"dash": "dash", "width": 1})

    # horizontal line for DFT final energy
    anno = {"text": "DFT final energy", "yanchor": "top"}
    fig.add_hline(
        y=dft_energy,
        line=dict(dash="dot", width=1, color=line_color),
        annotation=anno,
    )

    return fig


app = Dash(prevent_initial_callbacks=True, assets_folder=SETTINGS.ASSETS_PATH)

app_div = html.Div(
    [
        html.H1("Structure Relaxation Trajectory", style={"fontSize": "2em"}),
        html.P("Drag slider to see structure at different relaxation steps."),
    ],
    style=dict(margin="auto", textAlign="center", maxWidth="1200px", padding="2em"),
)

for name, df, traj in (
    ("CHGNet", df_chgnet, chgnet_traj),
    ("M3GNet", df_m3gnet, m3gnet_traj),
):
    step_size = max(1, len(df) // 20)  # ensure slider has max 20 steps
    slider = dcc.Slider(
        id=f"{name}-slider", min=0, max=len(df) - 1, step=step_size, updatemode="drag"
    )

    # don't create layout twice, causes duplicate
    # id error when restarting app in Jupyter notebook
    if name not in struct_layouts:
        struct_comps[name] = ctc.StructureMoleculeComponent(
            id=f"{name}-structure",
            struct_or_mol=init_struct,
        )
        struct_layouts[name] = struct_comps[name].layout()

    spg = init_struct.get_space_group_info()
    title = f"{material_id} - Spacegroup = {spg}"
    graph = dcc.Graph(
        id=f"{name}-fig",
        figure=plot_energy_and_forces(df, 0, e_col, force_col, title),
        style={"maxWidth": "50%"},
    )

    app_div.children += (
        html.H2(name, style=dict(margin="2em 0 1em", fontSize="1.5em")),
        slider,
        html.Div(
            [struct_layouts[name], graph],
            style=dict(display="flex", gap="2em", placeContent="center"),
        ),
    )

    def update_factory(
        trajectory: TrajectoryObserver,
        df_traj: pd.DataFrame,
        material_id: str,
        init_struct: Structure,
        e_col: str,
        force_col: str,
    ) -> Callable[[int], tuple[Structure, go.Figure]]:
        """Factory function to create update callback for each structure component.

        Args:
            trajectory (TrajectoryObserver): Contains relaxation trajectory.
            df_traj (pd.DataFrame): DataFrame with energy and force data.
            material_id (str): Material ID.
            init_struct (Structure): Initial structure.
            e_col (str): Column name for energy data.
            force_col (str): Column name for force data.

        Returns:
            Callable[[int], tuple[Structure, go.Figure]]: Update callback.
        """

        def update_structure(step: int) -> tuple[Structure, go.Figure]:
            """Update the structure displayed in the StructureMoleculeComponent and the
            dashed vertical line in the figure when the slider is moved.
            """
            lattice = trajectory.cells[step]
            coords = trajectory.atom_positions[step]
            init_struct.lattice = lattice
            if len(init_struct) != len(coords):
                raise ValueError(f"{len(init_struct)} != {len(coords)}")
            for idx, site in enumerate(init_struct):
                site.coords = coords[idx]

            spg = init_struct.get_space_group_info()
            title = f"{material_id} - Spacegroup = {spg}"

            fig = plot_energy_and_forces(df_traj, step, e_col, force_col, title)
            return init_struct, fig

        return update_structure

    app.callback(
        Output(struct_comps[name].id(), "data"),
        Output(graph, "figure"),
        Input(slider, "value"),
    )(update_factory(traj, df, material_id, init_struct, e_col, force_col))


app.layout = app_div
ctc.register_crystal_toolkit(app=app, layout=app.layout)
app.run(use_reloader=not IS_IPYTHON)
