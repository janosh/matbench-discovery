# %%
import gzip
import os
import pickle
import urllib.request
import warnings
from datetime import datetime
from glob import glob

import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from pymatgen.core import Structure
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility as MP2020Compat,
)
from pymatgen.entries.compatibility import (
    MaterialsProjectCompatibility as MPLegacyCompat,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz import density_scatter
from tqdm import tqdm

from mb_discovery import ROOT
from mb_discovery.energy import get_e_form_per_atom

try:
    import gdown
except ImportError:
    print(
        "gdown not installed. Needed for downloading WBM initial + relaxed structures "
        "from Google Drive."
    )

"""
Dataset generated with DFT and published in Jan 2021 as
"Predicting stable crystalline compounds using chemical similarity"
https://nature.com/articles/s41524-020-00481-6
"""


module_dir = os.path.dirname(__file__)
today = f"{datetime.now():%Y-%m-%d}"


# %% links to google drive files received via email from 1st author Hai-Chen Wang
# on 2021-06-15 containing initial and relaxed structures
google_drive_ids = {
    1: "1ZUgtYwrfZn_P8bULWRtTXepyAxHVxS5C",
    2: "1-3uu2AcARJxH7GReteGVASZTuttFGiW_",
    3: "1hc5BvDiFfTu_tc5F8m7ONSw2OgL9vN6o",
    4: "1aMYxG5YJUgMHpbWmHpzL4hRfmP26UQqh",
    5: "17kQt2r78ReWle4PhEIOXG7w7BFdezGM1",
}


# %%
for step, file_id in google_drive_ids.items():
    file_path = f"{module_dir}/raw/wbm-structures-step-{step}.json.bz2"

    if os.path.exists(file_path):
        print(f"{step=}: file exists, skipping {file_path=}")
        continue

    print(f"{step=}")
    gdown.download(f"https://drive.google.com/u/0/uc?id={file_id}", file_path)


# %%
summary_id_file = "1639IFUG7poaDE2uB6aISUOi65ooBwCIg"
summary_path = f"{module_dir}/raw/wbm-summary.txt"
gdown.download(f"https://drive.google.com/u/0/uc?id={summary_id_file}", summary_path)


# %%
json_paths = sorted(glob(f"{module_dir}/raw/wbm-structures-step-*.json.bz2"))
step_lens = (61848, 52800, 79205, 40328, 23308)
# step 3 has 79,211 structures but only 79,205 ComputedStructureEntries
# the 6 extra structures have missing energy, volume, etc. in the summary file

assert len(json_paths) == len(step_lens), "Mismatch in WBM steps and JSON files"
wbm_struct_json_checksums = (
    -7815922250032563359,
    -86268461085685423,
    -7707371069320539066,
    -3579196048285845088,
    -248039116266365352,
)


dfs_wbm_structs = {}
for json_path in json_paths:
    step = int(json_path.split(".json.bz2")[0][-1])
    assert step in range(1, 6)

    if step in dfs_wbm_structs:
        print(f"{step=} {json_path=} already loaded.")
        continue

    print(f"{step=}")
    df = pd.read_json(json_path).T

    # we hash index only for speed
    # could use joblib.hash(df) to hash whole df but it's slow
    checksum = pd.util.hash_pandas_object(df.index).sum()
    assert checksum == wbm_struct_json_checksums[step - 1], "bad JSON file checksum"

    if step == 3:
        # step 3 has 6 extra entries, see comment above
        additional_ids = (70802, 70803, 70825, 70826, 70828, 70829)
        df = df.drop(index=[f"step_3_{id}" for id in additional_ids])

    step_len = step_lens[step - 1]
    assert len(df) == step_len, f"bad len for {step=}: {len(df)} != {step_len}"
    dfs_wbm_structs[step] = df


# %%
df_all_steps = pd.concat(dfs_wbm_structs.values())

assert len(df_all_steps) == sum(step_lens)


def increment_wbm_material_id(wbm_id: str) -> str:
    """Maps step_1_0, step_1_1, ... onto wbm-step-1-1, wbm-step-1-2, ..."""
    prefix, step_num, material_num = wbm_id.split("_")

    assert prefix == "step"
    assert step_num.isdigit() and material_num.isdigit()

    return f"wbm-step-{step_num}-{int(material_num) + 1}"


df_all_steps.index = df_all_steps.index.map(increment_wbm_material_id)
df_all_steps.index.name = "material_id"
assert df_all_steps.index[0] == "wbm-step-1-1"
assert df_all_steps.index[-1] == "wbm-step-5-23308"

df_all_steps["initial_structure"] = df_all_steps.pop("org")
df_all_steps["final_structure"] = df_all_steps.pop("opt")
assert df_all_steps.columns == ["initial_structure", "final_structure"]

# df_all_steps = pd.read_json(
#     f"{module_dir}/2022-10-18-wbm-init+final-structures.json.bz2"
# )


# %% download WBM ComputedStructureEntries from
# https://archive.materialscloud.org/record/2021.68
mat_cloud_url = "https://archive.materialscloud.org/record/file?record_id=840"

for filename in (
    "README.txt",
    # same file as summary_path above except this doesn't have ID column
    # "summary.txt.bz2",
    *(f"step_{step}.json.bz2" for step in range(1, 6)),
):
    file_path = f"{module_dir}/raw/wbm-cse-{filename.lower().replace('_', '-')}"
    if os.path.exists(file_path):
        continue

    urllib.request.urlretrieve(f"{mat_cloud_url}&{filename=}", file_path)


# %%
cse_step_paths = sorted(glob(f"{module_dir}/raw/wbm-cse-step-*.json.bz2"))
assert len(cse_step_paths) == 5


"""
There is a discrepancy of 6 entries between the files on Materials Cloud containing the
ComputedStructureEntries (CSE) and those on Google Drive containing initial+relaxed
structures. The CSE files contain 257,489 entries but the summary file on materials
cloud and the org/opt (initial/relaxed) structures both have 257,495 rows. The
additional structures seem to be "step_3_70802", "step_3_70803", "step_3_70825",
"step_3_70826", "step_3_70828", and "step_3_70829". By removing these entries, the
compositions between the Materials Cloud and Google Drive data align.
"""

dfs_wbm_cses = {}
for json_path in tqdm(cse_step_paths):
    step = int(json_path.split(".json.bz2")[0][-1])

    if step in dfs_wbm_cses:
        print(f"{json_path=} already loaded.")
        continue

    df = pd.read_json(json_path)

    step_len = step_lens[step - 1]
    dfs_wbm_cses[step] = df
    assert len(df) == step_len, f"{step=}: {len(df)} != {step_len}"


# %%
df_all_steps["computed_structure_entry"] = pd.concat(dfs_wbm_cses.values()).to_numpy()

for cse in tqdm(df_all_steps.computed_structure_entry):
    # needed to ensure MaterialsProjectCompatibility can process the entries
    cse["parameters"]["run_type"] = (
        "GGA+U" if cse["parameters"]["is_hubbard"] else "GGA"
    )


assert pd.Series(
    cse["parameters"]["run_type"] for cse in tqdm(df_all_steps.computed_structure_entry)
).value_counts().to_dict() == {"GGA": 248481, "GGA+U": 9008}


# %%
# get composition from CSEs

df_all_steps["composition_from_cse"] = [
    ComputedStructureEntry.from_dict(cse).composition
    for cse in tqdm(df_all_steps.computed_structure_entry)
]

df_all_steps["composition_from_relaxed_struct"] = [
    Structure.from_dict(struct).composition
    for struct in tqdm(df_all_steps.final_structure)
]


# all but 1 composition matches between CSE and relaxed structure
# mismatching ID: wbm-step-1-37977 which becomes equal on reduction:
# CSE Comp: Ag4 Bi4 O12
# relaxed Comp: Ag16 Bi16 O48


df_unmatched = df_all_steps.query(
    "composition_from_cse != composition_from_relaxed_struct"
)
assert len(df_unmatched) == 1
assert df_unmatched.index[0] == "wbm-step-1-37977"
assert (
    df_unmatched.iloc[0].composition_from_cse.reduced_composition
    == df_unmatched.iloc[0].composition_from_relaxed_struct.reduced_composition
)

df_all_steps.pop("composition_from_relaxed_struct")  # not needed anymore


# %%
# randomly sample structures and ensure they match between CSE and relaxed structure
n_samples = 1000
for row in tqdm(df_all_steps.sample(n_samples).itertuples(), total=n_samples):
    struct_final = Structure.from_dict(row.final_structure)
    struct_from_cse = Structure.from_dict(row.computed_structure_entry["structure"])
    assert struct_final.matches(struct_from_cse), f"structure mismatch for {row.Index=}"

    # and check initial and final compositions match
    struct_init = Structure.from_dict(row.initial_structure)
    assert (
        struct_init.composition == struct_final.composition
    ), f"composition mismatch for {row.Index=}"


# %%
for mat_id, cse in df_all_steps.computed_structure_entry.items():
    cse["entry_id"] = mat_id

assert not any(
    entry["entry_id"] is None for entry in df_all_steps.computed_structure_entry
)

df_all_steps["formula_from_cse"] = [
    x.formula for x in df_all_steps.pop("composition_from_cse")
]
df_all_steps.drop(columns=["final_structure", "cse"]).to_json(
    f"{module_dir}/{today}-wbm-cses+init-structs.json.bz2"
)

# df_all_steps = pd.read_json(f"{module_dir}/2022-10-19-wbm-cses+init-structs.json.bz2")


# %%
col_map = {
    "# comp": "formula",
    "nsites": "n_sites",
    "vol": "volume",
    "e": "uncorrected_energy",
    "e_form": "e_form_per_atom",
    "e_hull": "e_hull",
    "gap": "bandgap_pbe",
    "id": "material_id",
}
# WBM summary was shared twice, once on google drive, once on materials cloud
# download both and check for consistency
df_summary = pd.read_csv(
    f"{module_dir}/raw/wbm-summary.txt", sep="\t", names=col_map.values()
).set_index("material_id")

df_summary_bz2 = pd.read_csv(
    f"{mat_cloud_url}&filename=summary.txt.bz2", sep="\t"
).rename(columns=col_map)

# duplicate Ga3Ru2U3 step_3_28147 (1st one is wbm-step-2-18689) has 0 volume in
# df_summary_bz2 vs 155.41 in df_summary
query_str = "volume > 0 & formula != 'Ga3Ru2U3'"
pd.testing.assert_frame_equal(
    df_summary.reset_index(drop=True).query(query_str),
    df_summary_bz2.reset_index(drop=True).query(query_str),
)

# fix bad energy in which is 0 in df_summary, -63.68 in CSE
df_summary.at["wbm-step-2-18689", "energy"] = df_all_steps.loc[
    "wbm-step-2-18689"
].computed_structure_entry["energy"]


# make sure dropping materials with 0 volume removes exactly 6 materials, the same ones
# listed in additional_ids above
assert len(df_summary.query("volume > 0")) == len(df_all_steps)
assert all(
    df_summary.reset_index().query("volume == 0").index.values - sum(step_lens[:2])
    == additional_ids
)
df_summary = df_summary.query("volume > 0")
df_summary.index = df_all_steps.index


# %% scatter plot summary energies vs CSE energies
df_summary["uncorrected_energy_from_cse"] = [
    cse["energy"] for cse in tqdm(df_all_steps.computed_structure_entry)
]

# check CSE and summary energies are consistent, only exceeding 0.1 eV difference twice
e_diff_summary_vs_cse = (
    df_summary.uncorrected_energy - df_summary.uncorrected_energy_from_cse
)
assert (
    e_diff_summary_vs_cse.max() < 0.15 and sum(e_diff_summary_vs_cse > 0.1) == 2
), df_summary.query("energy - uncorrected_energy_from_cse > 0.1")

density_scatter(df_summary.uncorrected_energy, df_summary.uncorrected_energy_from_cse)


# %%
df_all_steps["cse"] = [
    ComputedStructureEntry.from_dict(x)
    for x in tqdm(df_all_steps.computed_structure_entry)
]


# raw WBM ComputedStructureEntries have no energy corrections applied:
assert all(cse.uncorrected_energy == cse.energy for cse in df_all_steps.cse)
# summary and CSE n_sites match
assert all(df_summary.n_sites == len(cse.structure) for cse in df_all_steps.cse)


mp_compat = MP2020Compat()
compat_out = mp_compat.process_entries(df_all_steps.cse, clean=True, verbose=True)
assert len(compat_out) == len(df_all_steps) == len(df_summary)

n_corrected = sum(cse.uncorrected_energy != cse.energy for cse in df_all_steps.cse)
if isinstance(mp_compat, MPLegacyCompat):
    assert n_corrected == 39595, f"{n_corrected=}"
if isinstance(mp_compat, MP2020Compat):
    assert n_corrected == 100931, f"{n_corrected=}"

corr_label = "mp2020" if isinstance(mp_compat, MP2020Compat) else "legacy"
df_summary[f"e_corrections_{corr_label}"] = [
    cse.energy - cse.uncorrected_energy for cse in df_all_steps.cse
]


assert df_summary.e_corrections_mp2020.mean().round(4) == -0.9979
assert df_summary.e_corrections_legacy.mean().round(4) == -0.0643
assert (df_summary.filter(like="corrections").abs() > 1e-4).sum().to_dict() == {
    "e_corrections_mp2020": 100931,
    "e_corrections_legacy": 39595,
}


# %%
with gzip.open(
    f"{ROOT}/mb_discovery/energy/2022-10-13-rhys/ppd-mp.pkl.gz", "rb"
) as zip_file:
    ppd_rhys: PatchedPhaseDiagram = pickle.load(zip_file)


with gzip.open(f"{ROOT}/data/2022-09-18-ppd-mp.pkl.gz", "rb") as zip_file:
    ppd_mp = pickle.load(zip_file)


# %%
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

# takes ~20 min at 200 it/s for 250k entries in WBM
for entry in tqdm(df_all_steps.cse):
    assert entry.entry_id.startswith("wbm-step-")
    corr_label = "mp2020_" if isinstance(mp_compat, MP2020Compat) else "legacy_"
    # corr_label = "un"
    at_idx = entry.entry_id, f"e_above_hull_{corr_label}corrected_ppd_mp"

    if at_idx not in df_summary or pd.isna(df_summary.at[at_idx]):
        # use entry.(uncorrected_)energy_per_atom
        e_above_hull = (
            entry.corrected_energy_per_atom
            - ppd_mp.get_hull_energy_per_atom(entry.composition)
        )
        df_summary.at[at_idx] = e_above_hull


# %% calculate formation energies from CSEs wrt MP elemental reference energies
df_summary["e_form_per_atom_uncorrected"] = [
    get_e_form_per_atom(dict(composition=row.formula, energy=row.uncorrected_energy))
    for row in tqdm(df_summary.itertuples(), total=len(df_summary))
]

density_scatter(
    df=df_summary,
    x="e_form_per_atom_uncorrected",
    y="e_form_per_atom_mp2020_corrected",
)


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")


df_hull[df_all_steps.filter(like="e_above_hull").columns] = df_all_steps.filter(
    like="e_above_hull"
)

density_scatter(
    df=df_hull.query("e_above_hull_legacy_corrected < 3"),
    x="e_above_hull_mp",
    y="e_above_hull_legacy_corrected",
    xlabel=r"E$_\mathrm{above hull}$ from Rhys (legacy corrected)",
    ylabel=r"E$_\mathrm{above hull}$ from legacy corrected",
)


df_hull.query("abs(e_above_hull_mp2020_corrected - e_above_hull_mp) > 0.1")


# %%
df_summary.round(6).to_csv(f"{module_dir}/{today}-wbm-summary.csv")
