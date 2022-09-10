# %%
from __future__ import annotations

import gzip
import io
import pickle
import urllib.request
from datetime import datetime
from glob import glob

import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram, PDEntry
from pymatgen.core import Structure
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.plot_scripts.plot_funcs import (
    hist_classified_stable_as_func_of_hull_dist,
)

today = f"{datetime.now():%Y-%m-%d}"


# %%
task_type = "RS3RE"
date = "2022-08-19"
glob_pattern = f"{date}-m3gnet-wbm-relax-{task_type}/*.json.gz"
file_paths = sorted(glob(f"{ROOT}/data/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
# 2022-08-16 tried multiprocessing.Pool() to load files in parallel but was somehow
# slower than serial loading
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    try:
        # keep whole dataframe in memory
        df = pd.read_json(file_path)
        df.index = df.index.str.replace("_", "-")
        df.index.name = "material_id"
        col_map = dict(
            final_structure="m3gnet_structure", trajectory="m3gnet_trajectory"
        )
        df = df.rename(columns=col_map)
        df.reset_index().to_json(file_path)
        df["m3gnet_energy"] = df.m3gnet_trajectory.map(lambda x: x["energies"][-1][0])
        df["m3gnet_structure"] = df.m3gnet_structure.map(Structure.from_dict)
        df["formula"] = df.m3gnet_structure.map(lambda x: x.formula)
        df["volume"] = df.m3gnet_structure.map(lambda x: x.volume)
        df["n_sites"] = df.m3gnet_structure.map(len)
        dfs[file_path] = df.drop(columns=["m3gnet_trajectory"])
    except (ValueError, FileNotFoundError):
        # pandas v1.5+ correctly raises FileNotFoundError, below raises ValueError
        continue


# %%
df_m3gnet = pd.concat(dfs.values())
if any(df_m3gnet.index.str.contains("_")):
    df_m3gnet.index = df_m3gnet.index.str.replace("_", "-")


# %%
# 2022-01-25-ppd-mp+wbm.pkl.gz (235 MB)
ppd_pickle_url = "https://figshare.com/files/36669624"
zipped_file = urllib.request.urlopen(ppd_pickle_url)

ppd_mp_wbm: PatchedPhaseDiagram = pickle.load(
    io.BytesIO(gzip.decompress(zipped_file.read()))
)


pd_entries_m3gnet = [
    PDEntry(row.m3gnet_structure.composition, row.m3gnet_energy)
    for row in df_m3gnet.itertuples()
]
df_m3gnet["e_form_m3gnet_from_ppd"] = [
    ppd_mp_wbm.get_form_energy_per_atom(x) for x in pd_entries_m3gnet
]


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_m3gnet["e_above_mp_hull"] = df_hull.e_above_mp_hull


df_wbm = pd.read_csv(  # download wbm-steps-summary.csv (23.31 MB)
    "https://figshare.com/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")

df_m3gnet["e_form_wbm"] = df_wbm.e_form
df_m3gnet["wbm_energy"] = df_wbm.energy

pd_entries_wbm = [
    PDEntry(row.m3gnet_structure.composition, row.wbm_energy)
    for row in df_m3gnet.itertuples()
]
df_m3gnet["e_form_ppd_2022_01_25"] = [
    ppd_mp_wbm.get_form_energy_per_atom(x) for x in pd_entries_wbm
]


# %%
df_m3gnet.hist(bins=200, figsize=(18, 12))
df_m3gnet.isna().sum()


# %%
out_path = f"{ROOT}/data/{today}-m3gnet-wbm-relax-{task_type}.json.gz"
df_m3gnet.reset_index().to_json(out_path, default_handler=as_dict_handler)

out_path = f"{ROOT}/data/2022-08-16-m3gnet-wbm-relax-results-IS2RE.json.gz"
df_m3gnet = pd.read_json(out_path)


# %%
df_m3gnet["e_above_hull_pred"] = (  # TODO fix this incorrect e_above_hull_pred
    df_m3gnet["e_form_m3gnet_from_ppd"] - df_m3gnet["e_above_mp_hull"]
)

ax_hull_dist_hist = hist_classified_stable_as_func_of_hull_dist(
    e_above_hull_pred=df_m3gnet.e_above_hull_pred,
    e_above_hull_true=df_m3gnet.e_above_mp_hull,
)

# ax_hull_dist_hist.figure.savefig(f"{ROOT}/plots/{today}-m3gnet-wbm-hull-dist-hist.pdf")
