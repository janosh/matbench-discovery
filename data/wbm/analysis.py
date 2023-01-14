# %%
import os

import pandas as pd
from pymatviz import count_elements, ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import ROOT, today

module_dir = os.path.dirname(__file__)

"""
Compare MP and WBM elemental prevalence. Starting with WBM, MP below.
"""


# %%
df_summary = pd.read_csv(f"{module_dir}/2022-10-19-wbm-summary.csv").set_index(
    "material_id"
)
elem_counts = count_elements(df_summary.formula).astype(int)

elem_counts.to_json(
    f"{ROOT}/site/src/routes/about-the-test-set/{today}-wbm-element-counts.json"
)


# %%
fig = ptable_heatmap_plotly(
    elem_counts,
    log=True,
    colorscale="YlGnBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=elem_counts,
    font_size="1vw",
)

title = "WBM Elements"
fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
fig.show()


# %%
fig.write_image(f"{module_dir}/{today}-wbm-elements.svg", width=1000, height=500)
save_fig(fig, f"{module_dir}/{today}-wbm-elements.svelte")


# %% load MP training set
df = pd.read_json(f"{module_dir}/../mp/2022-08-13-mp-energies.json.gz")
elem_counts = count_elements(df.formula_pretty).astype(int)

elem_counts.to_json(
    f"{ROOT}/site/src/routes/about-the-test-set/{today}-mp-element-counts.json"
)
elem_counts.describe()


# %%
fig = ptable_heatmap_plotly(
    elem_counts[elem_counts > 1],
    log=True,
    colorscale="YlGnBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=elem_counts,
    font_size="1vw",
)

title = "MP Elements"
fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
fig.show()


# %%
fig.write_image(f"{module_dir}/{today}-mp-elements.svg", width=1000, height=500)
save_fig(fig, f"{module_dir}/{today}-mp-elements.svelte")
