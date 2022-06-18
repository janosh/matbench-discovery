# %%
import pandas as pd
import plotly.express as px
import plotly.io as pio
from diel_frontier.utils import decode_df_from_json

from ml_stability import ROOT


__author__ = "Janosh Riebesell"
__date__ = "2022-06-18"


pio.templates.default = "plotly_white"

pio.templates["plotly_white"]["layout"]["margin"] = dict(l=20, r=20, t=20, b=20)


# %%
df = pd.read_json(f"{ROOT}/data/wbm_cleaned.json.gz", orient="split")


# %%
df = decode_df_from_json(f"{ROOT}/data/mp.tar.bz2", orient="split")


# %%
df = decode_df_from_json(f"{ROOT}/data/wbm_cleaned.json.gz")


# %%
df_lattice = pd.json_normalize(df.cse.map(lambda x: x["structure"]["lattice"]))
df = df.join(df_lattice.add_prefix("final_lattice_"))

df_lattice = pd.json_normalize(df.initial_structure.map(lambda x: x["lattice"]))
df = df.join(df_lattice.add_prefix("initial_lattice_"))

df.isna().sum()
df.query("initial_lattice_matrix.isna()")


# %%
px.histogram(
    df.filter(like="volume"),
    nbins=1000,
    barmode="overlay",
    opacity=0.5,
    range_x=[0, 500],
)


# %% histogram of alpha lattice angles (similar results for beta and gamma)
fig = px.histogram(df.filter(like="alpha"), nbins=1000, barmode="overlay", log_y=True)
fig.write_image("plots/alpha_lattice_angles.png", scale=2)


# %%
px.histogram(
    df.filter(regex="_c$"), nbins=1000, log_y=True, barmode="overlay", opacity=0.5
)
