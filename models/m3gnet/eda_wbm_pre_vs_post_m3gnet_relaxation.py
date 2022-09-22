# %%
from datetime import datetime

import pandas as pd
import plotly.express as px
import plotly.io as pio
from pymatgen.core import Structure
from pymatgen.util.coord import pbc_diff
from pymatviz.utils import add_identity_line
from sklearn.metrics import r2_score

from mb_discovery import ROOT

__author__ = "Janosh Riebesell"
__date__ = "2022-06-18"


pio.templates.default = "plotly_white"

today = f"{datetime.now():%Y-%m-%d}"


# %%
df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
).set_index("material_id")


# %%
df_m3gnet_is2re = pd.read_json(
    f"{ROOT}/models/m3gnet/2022-08-16-m3gnet-wbm-relax-results-IS2RE.json.gz"
).set_index("material_id")
df_m3gnet_rs2re = pd.read_json(
    f"{ROOT}/models/m3gnet/2022-08-19-m3gnet-wbm-relax-results-RS2RE.json.gz"
).set_index("material_id")


# %% spread M3GNet post-pseudo-relaxation lattice params into separate columns
df_m3gnet_is2re["final_energy"] = df_m3gnet_is2re.trajectory.map(
    lambda x: x["energies"][-1][0]
)

df_m3gnet_lattice = pd.json_normalize(
    df_m3gnet_is2re.initial_structure.map(lambda x: x["lattice"])
).add_prefix("m3gnet_")
df_m3gnet_is2re[df_m3gnet_lattice.columns] = df_m3gnet_lattice.to_numpy()
df_m3gnet_is2re


# %% spread WBM initial and final lattice params into separate columns
df_m3gnet_is2re["final_wbm_structure"] = df_wbm.cse.map(lambda x: x["structure"])
df_wbm_final_lattice = pd.json_normalize(
    df_m3gnet_is2re.final_wbm_structure.map(lambda x: x["lattice"])
).add_prefix("final_wbm_")
df_m3gnet_is2re[df_wbm_final_lattice.columns] = df_wbm_final_lattice.to_numpy()


df_m3gnet_is2re["initial_wbm_structure"] = df_wbm.initial_structure
df_wbm_initial_lattice = pd.json_normalize(
    df_m3gnet_is2re.initial_structure.map(lambda x: x["lattice"])
).add_prefix("initial_wbm_")
df_m3gnet_is2re[df_wbm_initial_lattice.columns] = df_wbm_initial_lattice.to_numpy()


# %%
df_wbm_final_lattice = pd.json_normalize(
    df_wbm.cse.map(lambda x: x["structure"]["lattice"])
).add_prefix("final_wbm_")
df_wbm = df_wbm.join(df_wbm_final_lattice)

df_wbm_initial_lattice = pd.json_normalize(
    df_wbm.initial_structure.map(lambda x: x["lattice"])
).add_prefix("initial_wbm_")
df_wbm[df_wbm_initial_lattice.columns] = df_wbm_initial_lattice.to_numpy()

print(f"{df_wbm.isna().sum()=}")

df_wbm.query("initial_wbm_matrix.isna()")


# %%
px.histogram(
    df_m3gnet_is2re.filter(like="volume"),
    nbins=500,
    barmode="overlay",
    opacity=0.5,
    range_x=[0, 500],
)


# %%
fig = px.scatter(
    df_m3gnet_is2re.round(1),
    x="final_wbm_volume",
    y=["initial_wbm_volume", "m3gnet_volume"],
    hover_data=[df_m3gnet_is2re.index],
)
add_identity_line(fig)
fig.update_layout(
    title="Slightly tighter spread of M3GNet-relaxed vs initial WBM volumes"
)
fig.show()


# %% histogram of alpha lattice angles (similar results for beta and gamma)
fig = px.histogram(
    df_m3gnet_is2re.filter(like="alpha"), nbins=1000, barmode="overlay", log_y=True
)
fig.show()


# %%
px.histogram(
    df_m3gnet_is2re.filter(regex="_c$"),
    nbins=1000,
    log_y=True,
    barmode="overlay",
    opacity=0.5,
)


# %% compute mean absolute PBC difference between initial and final fractional
# coordinates of crystal sites
df_m3gnet_is2re["m3gnet_structure"] = df_m3gnet_is2re.m3gnet_structure.map(
    Structure.from_dict
)
df_m3gnet_is2re["initial_wbm_structure"] = df_m3gnet_is2re.initial_wbm_structure.map(
    Structure.from_dict
)
df_m3gnet_is2re["final_wbm_structure"] = df_m3gnet_is2re.final_wbm_structure.map(
    Structure.from_dict
)


df_m3gnet_is2re["m3gnet_pbc_diffs"] = [
    abs(
        pbc_diff(
            row.initial_wbm_structure.frac_coords,
            row.m3gnet_structure.frac_coords,
        )
    ).mean()
    for row in df_m3gnet_is2re.itertuples()
]


df_m3gnet_is2re["wbm_pbc_diffs"] = [
    abs(
        pbc_diff(
            row.initial_wbm_structure.frac_coords,
            row.final_wbm_structure.frac_coords,
        )
    ).mean()
    for row in df_m3gnet_is2re.itertuples()
]

df_m3gnet_is2re["m3gnet_to_final_wbm_pbc_diffs"] = [
    abs(
        pbc_diff(
            row.m3gnet_structure.frac_coords,
            row.final_wbm_structure.frac_coords,
        )
    ).mean()
    for row in df_m3gnet_is2re.itertuples()
]


print(
    "mean PBC difference of fractional coordinates before vs after relaxation with WBM "
    "and M3GNet"
)

wbm_pbc_diffs_mean = df_m3gnet_is2re.wbm_pbc_diffs.mean()
print(f"{wbm_pbc_diffs_mean = :.3}")

m3gnet_pbc_diffs_mean = df_m3gnet_is2re.m3gnet_pbc_diffs.mean()
print(f"{m3gnet_pbc_diffs_mean = :.3}")

m3gnet_to_final_wbm_pbc_diffs_mean = (
    df_m3gnet_is2re.m3gnet_to_final_wbm_pbc_diffs.mean()
)
print(f"{m3gnet_to_final_wbm_pbc_diffs_mean = :.3}")

print(f"{wbm_pbc_diffs_mean / m3gnet_pbc_diffs_mean = :.3}")


# %%
# plt_fig = df_m3gnet_is2re.plot.scatter(
#     x="e_m3gnet_per_atom_rs2re", y="e_m3gnet_per_atom_is2re"
# )
# df_m3gnet_is2re.filter(like="m3gnet_energy").hist(bins=100)

df_m3gnet_is2re["m3gnet_energy_rs2re"] = df_m3gnet_rs2re.m3gnet_energy

for task_type in ["is2re", "rs2re"]:
    energy_per_atom = (
        df_m3gnet_is2re[f"m3gnet_energy_{task_type}"] / df_m3gnet_is2re.n_sites
    )

    df_m3gnet_is2re[f"e_m3gnet_per_atom_{task_type}"] = energy_per_atom

fig = px.scatter(
    df_m3gnet_is2re,
    x="e_m3gnet_per_atom_rs2re",
    y="e_m3gnet_per_atom_is2re",
    render_mode="webgl",
)
add_identity_line(fig)

len_overlap = df_m3gnet_is2re.filter(like="e_m3gnet_per_atom_").dropna().shape[0]
x_vals, y_vals = df_m3gnet_is2re.filter(like="e_m3gnet_per_atom_").dropna().values.T

MAE = abs(x_vals - y_vals).mean()
R2 = r2_score(x_vals, y_vals)

title = f"data size = {len_overlap:,} \t {MAE = :.2} \t {R2 = :.4}"
fig.update_layout(title=dict(text=title, x=0.5))

# 250k scatter points require exporting to PNG, interactive version freezes the
# notebook server
fig.show(renderer="png", scale=2)
fig.write_image(
    f"{ROOT}/figures/{today}-m3gnet-energy-per-atom-scatter-is2re-vs-rs2re.png", scale=2
)


# %% write df back to compressed JSON
# filter out columns containing 'rs2re'
# df_m3gnet_is2re.reset_index().filter(regex="^((?!rs2re).)*$").to_json(
#     f"{ROOT}/models/m3gnet/2022-08-16-m3gnet-wbm-relax-results-IS2RE-2.json.gz"
# ).set_index("material_id")
