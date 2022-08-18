# %%
from __future__ import annotations

import gzip
import io
import pickle
from datetime import datetime
from glob import glob
from urllib.request import urlopen

import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram, PDEntry
from pymatgen.core import Structure
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.plot_scripts.plot_funcs import (
    hist_classify_stable_as_func_of_hull_dist,
)


today = f"{datetime.now():%Y-%m-%d}"


# %%
glob_pattern = "2022-08-16-m3gnet-wbm-relax-results/*.json.gz"
file_paths = glob(f"{ROOT}/data/{glob_pattern}")
print(f"Found {len(file_paths):,} files for {glob_pattern = }")


dfs: dict[str, pd.DataFrame] = {}
# 2022-08-16 tried multiprocessing.Pool() to load files in parallel but was somehow
# slower than serial loading
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    try:
        dfs[file_path] = pd.read_json(file_path)
    except (ValueError, FileNotFoundError):
        # pandas v1.5+ correctly raises FileNotFoundError, below raises ValueError
        continue


# %%
df_m3gnet = pd.concat(dfs.values())
df_m3gnet.index.name = "material_id"
if any(df_m3gnet.index.str.contains("_")):
    df_m3gnet.index = df_m3gnet.index.str.replace("_", "-")

df_m3gnet = df_m3gnet.rename(
    columns=dict(final_structure="m3gnet_structure", trajectory="m3gnet_trajectory")
)

df_m3gnet["m3gnet_energy"] = df_m3gnet.trajectory.map(lambda x: x["energies"][-1][0])


# %%
# 2022-01-25-ppd-mp+wbm.pkl.gz (235 MB)
ppd_pickle_url = "https://figshare.com/ndownloader/files/36669624"
zipped_file = urlopen(ppd_pickle_url)

ppd_mp_wbm: PatchedPhaseDiagram = pickle.load(
    io.BytesIO(gzip.decompress(zipped_file.read()))
)


df_m3gnet["m3gnet_structure"] = df_m3gnet.m3gnet_structure.map(Structure.from_dict)
df_m3gnet["pd_entry"] = [
    PDEntry(row.m3gnet_structure.composition, row.m3gnet_energy)
    for row in df_m3gnet.itertuples()
]
df_m3gnet["e_form_m3gnet"] = df_m3gnet.pd_entry.map(ppd_mp_wbm.get_form_energy_per_atom)


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_m3gnet["e_above_mp_hull"] = df_hull.e_above_mp_hull


df_summary = pd.read_csv(f"{ROOT}/data/wbm-steps-summary.csv", comment="#").set_index(
    "material_id"
)

df_m3gnet["e_form_wbm"] = df_summary.e_form


# %%
df_m3gnet.hist(bins=200, figsize=(18, 12))
df_m3gnet.isna().sum()


# %%
out_path = f"{ROOT}/data/{today}-m3gnet-wbm-relax-results.json.gz"
df_m3gnet.drop(columns=["pd_entry"]).reset_index().to_json(
    out_path, default_handler=as_dict_handler
)


# %%
ax_hull_dist_hist = hist_classify_stable_as_func_of_hull_dist(
    formation_energy_targets=df_m3gnet.e_form_wbm,
    formation_energy_preds=df_m3gnet.e_form_m3gnet,
    e_above_hull_vals=df_m3gnet.e_above_mp_hull,
)

# ax_hull_dist_hist.figure.savefig(f"{ROOT}/plots/{today}-m3gnet-wbm-hull-dist-hist.pdf")
