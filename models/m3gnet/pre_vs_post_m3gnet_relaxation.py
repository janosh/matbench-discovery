"""Compare M3GNet-relaxed vs DFT-relaxed WBM lattice volumes and angles."""

# %%
import os

import pandas as pd
import plotly.express as px
import pymatviz as pmv
from pymatgen.core import Structure
from pymatgen.util.coord import pbc_diff
from pymatviz.enums import Key
from pymatviz.powerups import add_identity_line
from sklearn.metrics import r2_score

from matbench_discovery import ROOT, SITE_FIGS, plots
from matbench_discovery.data import DataFiles

__author__ = "Janosh Riebesell"
__date__ = "2022-06-18"

module_dir = os.path.dirname(__file__)
del plots  # https://github.com/PyCQA/pyflakes/issues/366


# %%
df_wbm = pd.read_json(DataFiles.wbm_cses_plus_init_structs.path).set_index(Key.mat_id)

df_summary = pd.read_csv(DataFiles.wbm_summary.path).set_index(Key.mat_id)


# %%
is2re_path = f"{ROOT}/models/m3gnet/2022-10-31-m3gnet-wbm-IS2RE.json.gz"
df_m3gnet_is2re = pd.read_json(is2re_path).set_index(Key.mat_id)

rs2re_path = f"{ROOT}/models/m3gnet/2022-08-19-m3gnet-wbm-RS2RE.json.gz"
df_m3gnet_rs2re = pd.read_json(rs2re_path).set_index(Key.mat_id)


# %%
df_wbm["m3gnet_volume"] = df_m3gnet_is2re.m3gnet_volume


# %% spread M3GNet post-pseudo-relaxation lattice params into separate columns
df_m3gnet_lattice = pd.json_normalize(
    df_m3gnet_is2re[Key.init_struct].map(lambda x: x["lattice"])
).add_prefix("m3gnet_")
df_m3gnet_is2re[df_m3gnet_lattice.columns] = df_m3gnet_lattice.to_numpy()

# df_m3gnet_is2re["m3gnet_energy"] = df_m3gnet_is2re.trajectory.map(
#     lambda x: x["energies"][-1][0]
# )


# %% spread WBM initial and final lattice params into separate columns
df_wbm_final_lattice = pd.json_normalize(
    df_wbm[Key.cse].map(lambda cse: cse["structure"]["lattice"])
).add_prefix("final_wbm_")
df_wbm["final_wbm_volume"] = df_wbm_final_lattice.final_wbm_volume.to_numpy()

df_wbm_initial_lattice = pd.json_normalize(
    df_wbm[Key.init_struct].map(lambda x: (x or {}).get("lattice"))
).add_prefix("initial_wbm_")
df_wbm["initial_wbm_volume"] = df_wbm_initial_lattice.initial_wbm_volume.to_numpy()


# 2 materials have no initial structure: wbm-5-23166, wbm-5-23294
print(f"{df_wbm.isna().sum()=}")
df_wbm.query("initial_wbm_volume.isna()").index.tolist()


# %% parity plot of M3GNet/initial volumes vs DFT-relaxed volumes
ax = pmv.density_scatter(
    df=df_wbm.query("m3gnet_volume < 2000"),
    x="final_wbm_volume",
    y="m3gnet_volume",
    cmap="Reds",
    alpha=0.5,
    stats=dict(loc="lower right", prefix="m3gnet to final (red)\n"),
)
pmv.density_scatter(
    df=df_wbm.query("m3gnet_volume < 2000"),
    x="final_wbm_volume",
    y="initial_wbm_volume",
    ax=ax,
    cmap="Blues",
    alpha=0.5,
    stats=dict(loc="upper left", prefix="init to final (blue)\n"),
)
ax.set(title="M3GNet-relaxed vs DFT-relaxed WBM volumes")
ax.set(xlabel="DFT-relaxed volume [Å³]")
ax.set(ylabel="M3GNet-relaxed / unrelaxed volume [Å³]")
pmv.save_fig(ax, f"{SITE_FIGS}/m3gnet-wbm-volume-scatter.webp", dpi=200)


# %% histogram of M3GNet-relaxed vs initial WBM volume residuals wrt DFT-relaxed volume
df_plot = df_wbm.query("m3gnet_volume < 300").filter(like="volume")
df_plot["m3gnet_vol_diff"] = df_plot.m3gnet_volume - df_plot.final_wbm_volume
df_plot["dft_vol_diff"] = df_plot.initial_wbm_volume - df_plot.final_wbm_volume
fig = px.histogram(
    df_plot.melt(
        value_vars=["m3gnet", "dft"], value_name="vol_diff", var_name="method"
    ),
    x="vol_diff",
    color="method",
    range_x=[-50, 50],
    barmode="overlay",
)
fig.show()
fig.write_image(f"{SITE_FIGS}/m3gnet-wbm-volume-diff-residual-hist.webp", scale=2)


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
print(f"{wbm_pbc_diffs_mean=:.3}")

m3gnet_pbc_diffs_mean = df_m3gnet_is2re.m3gnet_pbc_diffs.mean()
print(f"{m3gnet_pbc_diffs_mean=:.3}")

m3gnet_to_final_wbm_pbc_diffs_mean = (
    df_m3gnet_is2re.m3gnet_to_final_wbm_pbc_diffs.mean()
)
print(f"{m3gnet_to_final_wbm_pbc_diffs_mean=:.3}")

print(f"{wbm_pbc_diffs_mean / m3gnet_pbc_diffs_mean=:.3}")


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
x_vals, y_vals = df_m3gnet_is2re.filter(like="e_m3gnet_per_atom_").dropna().to_numpy().T

MAE = abs(x_vals - y_vals).mean()
R2 = r2_score(x_vals, y_vals)

title = f"data size = {len_overlap:,} \t {MAE = :.2} \t {R2 = :.4}"
fig.update_layout(title=dict(text=title, x=0.5))

# 250k scatter points require exporting to PNG, interactive version freezes the
# notebook server
fig.show(renderer="png", scale=2)
fig.write_image(
    f"{SITE_FIGS}/m3gnet-energy-per-atom-parity-is2re-vs-rs2re.webp", scale=2
)


# %% write df back to compressed JSON
# filter out columns containing 'rs2re'
# df_m3gnet_is2re.reset_index().filter(regex="^((?!rs2re).)*$").to_json(
#     f"{module_dir}/2022-10-31-m3gnet-wbm-IS2RE-2.json.gz"
# ).set_index(Key.mat_id)
