"""Compare CHGNet long vs short relaxations."""


# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.core import Structure
from pymatviz import density_scatter, plot_structure_2d, ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import PDF_FIGS
from matbench_discovery import plots as plots
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.preds import PRED_FILES

__author__ = "Janosh Riebesell"
__date__ = "2023-03-06"

module_dir = os.path.dirname(__file__)
id_col = "material_id"


# %%
df_chgnet = pd.read_csv(PRED_FILES.CHGNet)
df_chgnet = df_chgnet.set_index(id_col).add_suffix("_2000")
df_chgnet_500 = pd.read_csv(PRED_FILES.CHGNet.replace("-06", "-04"))
df_chgnet_500 = df_chgnet_500.set_index(id_col).add_suffix("_500")
df_chgnet[list(df_chgnet_500)] = df_chgnet_500
df_chgnet["formula"] = df_wbm.formula

e_form_2000 = "e_form_per_atom_chgnet_2000"
e_form_500 = "e_form_per_atom_chgnet_500"

min_e_diff = 0.1
# structures with smaller energy after longer relaxation need many steps
df_long = df_chgnet.query(f"{e_form_2000} - {e_form_500} < -{min_e_diff}")
# structures with larger energy after longer relaxation are problematic
df_bad = df_chgnet.query(f"{e_form_2000} - {e_form_500} > {min_e_diff}")
# both combined
df_diff = df_chgnet.query(f"abs({e_form_2000} - {e_form_500}) > {min_e_diff}")

assert len(df_long) + len(df_bad) == len(df_diff)


# %%
density_scatter(df=df_chgnet, x=e_form_500, y=e_form_2000)


# %%
df_diff.reset_index().plot.scatter(
    x=e_form_500,
    y=e_form_2000,
    hover_name=id_col,
    hover_data=["formula"],
    backend="plotly",
    title=f"{len(df_diff)} structures have > {min_e_diff} eV/atom energy diff after "
    "longer relaxation",
)


# %%
fig = ptable_heatmap_plotly(df_bad.formula)
title = "structures with larger error<br>after longer relaxation"
fig.layout.title.update(text=f"{len(df_diff)} {title}", x=0.4, y=0.9)
fig.show()


# %%
df_cse = pd.read_json(DATA_FILES.wbm_cses_plus_init_structs).set_index(id_col)
df_cse.loc[df_diff.index].reset_index().to_json(
    f"{module_dir}/wbm-chgnet-bad-relax.json.gz"
)


# %%
n_rows, n_cols = 3, 4
fig, axs = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 4 * n_rows))
n_struct = min(n_rows * n_cols, len(df_diff))
struct_col = "initial_structure"

fig.suptitle(f"{n_struct} {struct_col} {title}", fontsize=16, fontweight="bold", y=1.05)
for idx, (ax, row) in enumerate(
    zip(axs.flat, df_cse.loc[df_diff.index].itertuples()), 1
):
    struct = Structure.from_dict(getattr(row, struct_col))
    plot_structure_2d(struct, ax=ax)
    _, spg_num = struct.get_space_group_info()
    formula = struct.composition.reduced_formula
    id = row.Index
    ax.set_title(f"{idx}. {formula} (spg={spg_num})\n{id}", fontweight="bold")

save_fig(fig, f"{PDF_FIGS}/chgnet-bad-relax-structures.pdf")
