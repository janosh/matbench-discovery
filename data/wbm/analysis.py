# %%
import os

import pandas as pd
from pymatviz import count_elements, ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import df_wbm

module_dir = os.path.dirname(__file__)

"""
Compare MP and WBM elemental prevalence. Starting with WBM, MP below.
"""


# %%
wbm_elem_counts = count_elements(df_wbm.formula).astype(int)

# wbm_elem_counts.to_json(
#     f"{ROOT}/site/src/routes/about-the-test-set/{today}-wbm-element-counts.json"
# )


# %%
wbm_fig = ptable_heatmap_plotly(
    wbm_elem_counts.drop("Xe"),
    log=True,
    colorscale="RdBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=wbm_elem_counts,
)

title = "WBM Elements"
wbm_fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
wbm_fig.show()


# %%
wbm_fig.write_image(
    f"{module_dir}/figs/{today}-wbm-elements.svg", width=1000, height=500
)
save_fig(wbm_fig, f"{FIGS}/{today}-wbm-elements.svelte")


# %% load MP training set
df = pd.read_json(f"{module_dir}/../mp/2022-08-13-mp-energies.json.gz")
mp_elem_counts = count_elements(df.formula_pretty).astype(int)

# mp_elem_counts.to_json(
#     f"{ROOT}/site/src/routes/about-the-test-set/{today}-mp-element-counts.json"
# )
mp_elem_counts.describe()


# %%
mp_fig = ptable_heatmap_plotly(
    mp_elem_counts[mp_elem_counts > 1],
    log=True,
    colorscale="RdBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=mp_elem_counts,
)

title = "MP Elements"
mp_fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
mp_fig.show()


# %%
mp_fig.write_image(f"{module_dir}/figs/{today}-mp-elements.svg", width=1000, height=500)
# save_fig(mp_fig, f"{FIGS}/{today}-mp-elements.svelte")
