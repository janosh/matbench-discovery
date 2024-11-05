"""Compare CHGNet long vs short relaxations."""

# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
import pymatviz as pmv
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS
from matbench_discovery import plots as plots
from matbench_discovery.data import DataFiles, Model, df_wbm
from matbench_discovery.preds import df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-03-06"

module_dir = os.path.dirname(__file__)


# %%
df_chgnet = df_chgnet_v030 = pd.read_csv(Model.chgnet.path)
df_chgnet_v020 = pd.read_csv(
    f"{module_dir}/2023-03-06-chgnet-0.2.0-wbm-IS2RE.csv.gz", index_col=Key.mat_id
)
df_chgnet[Key.formula] = df_wbm[Key.formula]

e_form_2000 = "e_form_per_atom_chgnet_relax_steps_2000"
e_form_500 = "e_form_per_atom_chgnet_relax_steps_500"

min_e_diff = 0.1
# structures with smaller energy after longer relaxation need many steps
df_long = df_chgnet.query(f"{e_form_2000} - {e_form_500} < -{min_e_diff}")
# structures with larger energy after longer relaxation are problematic
df_bad = df_chgnet.query(f"{e_form_2000} - {e_form_500} > {min_e_diff}")
# both combined
df_diff = df_chgnet.query(f"abs({e_form_2000} - {e_form_500}) > {min_e_diff}")

if len(df_long) + len(df_bad) != len(df_diff):
    raise ValueError(f"{len(df_long)=} + {len(df_bad)=} != {len(df_diff)=}")


# %%
pmv.density_scatter(df=df_chgnet, x=e_form_500, y=e_form_2000)


# %%
df_diff.reset_index().plot.scatter(
    x=e_form_500,
    y=e_form_2000,
    hover_name=Key.mat_id,
    hover_data=[Key.formula],
    backend=pmv.utils.PLOTLY,
    title=f"{len(df_diff)} structures have > {min_e_diff} eV/atom energy diff after "
    "longer relaxation",
)


# %%
fig = pmv.ptable_heatmap_plotly(df_bad[Key.formula])
title = "structures with larger error<br>after longer relaxation"
fig.layout.title.update(text=f"{len(df_diff)} {title}", x=0.4, y=0.9)
fig.show()


# %%
df_cse = pd.read_json(DataFiles.wbm_cses_plus_init_structs.path).set_index(Key.mat_id)
df_cse.loc[df_diff.index].reset_index().to_json(
    f"{module_dir}/wbm-chgnet-bad-relax.json.gz"
)


# %%
n_rows, n_cols = 3, 4
fig, axs = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 4 * n_rows))
n_struct = min(n_rows * n_cols, len(df_diff))
struct_col = Key.init_struct

fig.suptitle(f"{n_struct} {struct_col} {title}", fontsize=16, fontweight="bold", y=1.05)
for idx, row in enumerate(df_cse.loc[df_diff.index].itertuples(), start=1):
    struct = Structure.from_dict(getattr(row, struct_col))
    ax = pmv.structure_2d(struct, ax=axs.flat[idx - 1])
    _, spg_num = struct.get_space_group_info()
    formula = struct.composition.reduced_formula
    ax.set_title(f"{idx}. {formula} (spg={spg_num})\n{row.Index}", fontweight="bold")

pmv.save_fig(fig, f"{PDF_FIGS}/chgnet-bad-relax-structures.pdf")


# %% ensure all CHGNet static predictions (direct energy without any structure
# relaxation) are higher in energy than the relaxed ones, i.e. that the optimizer is
# working correctly
pmv.density_scatter(df=df_preds, x="CHGNet", y="chgnet_no_relax")
