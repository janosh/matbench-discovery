# %%
import gzip
import os
import pickle
import urllib.request
import warnings
from datetime import datetime
from glob import glob

import pandas as pd
from aviary.wren.utils import get_aflow_label_from_spglib
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from pymatgen.core import Composition, Structure
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility as MP2020Compat,
)
from pymatgen.entries.compatibility import (
    MaterialsProjectCompatibility as MPLegacyCompat,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz import density_scatter
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.energy import get_e_form_per_atom

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
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


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
        print(f"{file_path} already exists, skipping")
        continue

    print(f"{step=}")
    gdown.download(f"https://drive.google.com/u/0/uc?id={file_id}", file_path)


# %%
summary_path = f"{module_dir}/raw/wbm-summary.txt"

if not os.path.exists(summary_path):
    summary_id_file = "1639IFUG7poaDE2uB6aISUOi65ooBwCIg"
    summary_url = f"https://drive.google.com/u/0/uc?id={summary_id_file}"
    gdown.download(summary_url, summary_path)


# %%
json_paths = sorted(glob(f"{module_dir}/raw/wbm-structures-step-*.json.bz2"))
step_lens = (61848, 52800, 79205, 40328, 23308)
# step 3 has 79,211 structures but only 79,205 ComputedStructureEntries
# i.e. 6 extra structures which have missing energy, volume, etc. in the summary file
bad_struct_ids = (70802, 70803, 70825, 70826, 70828, 70829)


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
        df = df.drop(index=[f"step_3_{id}" for id in bad_struct_ids])
        # re-index after dropping bad structures to get same indices as summary file
        # where IDs are consecutive, i.e. step_3_70801 is followed by step_3_70802,
        # not step_3_70804, etc.
        df.index = [f"step_3_{idx + 1}" for idx in range(len(df))]

    step_len = step_lens[step - 1]
    assert len(df) == step_len, f"bad len for {step=}: {len(df)} != {step_len}"
    dfs_wbm_structs[step] = df


# NOTE step 5 is missing 2 initial structures
assert dict(dfs_wbm_structs[5].isna().sum()) == {"opt": 0, "org": 2}
assert list(dfs_wbm_structs[5].query("org.isna()").index) == [
    "step_5_23165",
    "step_5_23293",
]


# %%
df_wbm = pd.concat(dfs_wbm_structs.values())

assert len(df_wbm) == sum(step_lens)


def increment_wbm_material_id(wbm_id: str) -> str:
    """Maps step_1_0, step_1_1, ... onto wbm-step-1-1, wbm-step-1-2, ..."""
    try:
        prefix, step_num, material_num = wbm_id.split("_")
    except ValueError:
        print(f"bad {wbm_id=}")
        return wbm_id

    msg = f"bad {wbm_id=}, {prefix=} {step_num=} {material_num=}"
    assert prefix == "step", msg
    assert step_num.isdigit(), msg
    assert material_num.isdigit(), msg

    return f"wbm-step-{step_num}-{int(material_num) + 1}"


df_wbm.index = df_wbm.index.map(increment_wbm_material_id)
df_wbm.index.name = "material_id"
assert df_wbm.index[0] == "wbm-step-1-1"
assert df_wbm.index[-1] == "wbm-step-5-23308"

df_wbm["initial_structure"] = df_wbm.pop("org")
df_wbm["final_structure"] = df_wbm.pop("opt")
assert list(df_wbm.columns) == ["initial_structure", "final_structure"]


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
        print(f"{file_path} already exists, skipping")
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
for json_path in cse_step_paths:
    step = int(json_path.split(".json.bz2")[0][-1])
    print(f"{step=}")
    assert step in range(1, 6)

    if step in dfs_wbm_cses:
        print(f"{json_path=} already loaded.")
        continue

    df = pd.read_json(json_path)

    step_len = step_lens[step - 1]
    dfs_wbm_cses[step] = df
    assert len(df) == step_len, f"{step=}: {len(df)} != {step_len}"


# %%
df_wbm["computed_structure_entry"] = pd.concat(dfs_wbm_cses.values()).to_numpy()

for mat_id, cse in df_wbm.computed_structure_entry.items():
    # needed to ensure MaterialsProjectCompatibility can process the entries
    cse["parameters"]["run_type"] = (
        "GGA+U" if cse["parameters"]["is_hubbard"] else "GGA"
    )
    cse["entry_id"] = mat_id
    assert cse["entry_id"].startswith("wbm-step-")

assert pd.Series(
    cse["parameters"]["run_type"] for cse in tqdm(df_wbm.computed_structure_entry)
).value_counts().to_dict() == {"GGA": 248481, "GGA+U": 9008}


# %% get composition from CSEs
df_wbm["composition_from_cse"] = [
    ComputedStructureEntry.from_dict(cse).composition
    for cse in tqdm(df_wbm.computed_structure_entry)
]

df_wbm["composition_from_final_struct"] = [
    Structure.from_dict(struct).composition for struct in tqdm(df_wbm.final_structure)
]

# all but 1 composition matches between CSE and final structure
# mismatching ID: wbm-step-1-37977 which becomes equal on reduction:
# CSE Comp: Ag4 Bi4 O12
# final structure Comp: Ag16 Bi16 O48
df_mismatch = df_wbm.query("composition_from_cse != composition_from_final_struct")
assert len(df_mismatch) == 1
assert df_mismatch.index[0] == "wbm-step-1-37977"
assert (
    df_mismatch.iloc[0].composition_from_cse.reduced_composition
    == df_mismatch.iloc[0].composition_from_final_struct.reduced_composition
)

df_wbm.pop("composition_from_final_struct")  # not needed anymore


# %% randomly sample structures and ensure they match between CSE and final structure
n_samples = 1000
for row in tqdm(df_wbm.sample(n_samples).itertuples(), total=n_samples):
    struct_final = Structure.from_dict(row.final_structure)
    struct_from_cse = Structure.from_dict(row.computed_structure_entry["structure"])
    assert struct_final.matches(struct_from_cse), f"structure mismatch for {row.Index=}"

    # and check initial and final compositions match
    struct_init = Structure.from_dict(row.initial_structure)
    assert (
        struct_init.composition == struct_final.composition
    ), f"composition mismatch for {row.Index=}"


# %%
df_wbm["formula_from_cse"] = [
    x.alphabetical_formula for x in df_wbm.pop("composition_from_cse")
]

for key, col_name in (
    ("cses", "computed_structure_entry"),
    ("init-structs", "initial_structure"),
):
    cols = ["initial_structure", "formula_from_cse", col_name]
    df_wbm[cols].reset_index().to_json(f"{module_dir}/{today}-wbm-{key}.json.bz2")


# %%
col_map = {
    "# comp": "formula",
    "nsites": "n_sites",
    "vol": "volume",
    "e": "uncorrected_energy",
    "e_form": "e_form_per_atom_wbm",
    "e_hull": "e_hull_wbm",
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


# make sure dropping materials with 0 volume removes exactly 6 materials, the same ones
# listed in bad_struct_ids above
assert len(df_summary.query("volume > 0")) == len(df_wbm)
assert all(
    df_summary.reset_index().query("volume == 0").index.values - sum(step_lens[:2])
    == bad_struct_ids
)
df_summary = df_summary.query("volume > 0")
df_summary.index = df_summary.index.map(increment_wbm_material_id)
assert sum(df_summary.index != df_wbm.index) == 0

# sort formulas alphabetically
df_summary["alph_formula"] = [
    Composition(x).alphabetical_formula for x in df_summary.formula
]
assert sum(df_summary.alph_formula != df_summary.formula) == 219_215
assert df_summary.alph_formula[3] == "Ag2 Au1 Hg1"
assert df_summary.formula[3] == "Ag2 Hg1 Au1"

df_summary["formula"] = df_summary.pop("alph_formula")


# %%
# check summary and CSE formulas agree
assert all(df_summary["formula"] == df_wbm.formula_from_cse)


# fix bad energy which is 0 in df_summary but a more realistic -63.68 in CSE
df_summary.at["wbm-step-2-18689", "uncorrected_energy"] = df_wbm.loc[
    "wbm-step-2-18689"
].computed_structure_entry["energy"]

# NOTE careful with ComputedEntries as object vs as dicts, the meaning of keys changes:
# cse.energy == cse.uncorrected_energy + cse.correction
# whereas
# cse.as_dict()["energy"] == cse.uncorrected_energy


# %% scatter plot summary energies vs CSE energies
df_summary["uncorrected_energy_from_cse"] = [
    cse["energy"] for cse in tqdm(df_wbm.computed_structure_entry)
]

# check CSE and summary energies are consistent, only exceeding 0.1 eV difference twice
diff_e_cse_e_summary = (
    df_summary.uncorrected_energy - df_summary.uncorrected_energy_from_cse
)
assert (
    diff_e_cse_e_summary.max() < 0.15 and sum(diff_e_cse_e_summary > 0.1) == 2
), df_summary.query("energy - uncorrected_energy_from_cse > 0.1")

density_scatter(df_summary.uncorrected_energy, df_summary.uncorrected_energy_from_cse)


# %%
# raw WBM ComputedStructureEntries have no energy corrections applied:
assert all(cse.uncorrected_energy == cse.energy for cse in df_wbm.cse)
# summary and CSE n_sites match
assert all(df_summary.n_sites == [len(cse.structure) for cse in df_wbm.cse])


mp_compat = MP2020Compat() if False else MPLegacyCompat()
compat_out = mp_compat.process_entries(df_wbm.cse, clean=True, verbose=True)

mp_compat.process_entry(cse)
assert len(compat_out) == len(df_wbm) == len(df_summary)

n_corrected = sum(cse.uncorrected_energy != cse.energy for cse in df_wbm.cse)
if isinstance(mp_compat, MPLegacyCompat):
    assert n_corrected == 39595, f"{n_corrected=}"
if isinstance(mp_compat, MP2020Compat):
    assert n_corrected == 100931, f"{n_corrected=}"

corr_label = "mp2020" if isinstance(mp_compat, MP2020Compat) else "legacy"
df_summary[f"e_correction_per_atom_{corr_label}"] = [
    cse.correction_per_atom for cse in df_wbm.cse
]

assert df_summary.e_correction_per_atom_mp2020.mean().round(4) == -0.1067
assert df_summary.e_correction_per_atom_legacy.mean().round(4) == -0.0643
assert (df_summary.filter(like="correction").abs() > 1e-4).sum().to_dict() == {
    "e_correction_per_atom_mp2020": 100931,
    "e_correction_per_atom_legacy": 39595,
}, "unexpected number of materials received non-zero corrections"

ax = density_scatter(
    df_summary.e_correction_per_atom_legacy,
    df_summary.e_correction_per_atom_mp2020,
    xlabel="legacy corrections (eV / atom)",
    ylabel="MP2020 corrections (eV / atom)",
)
# ax.figure.savefig(f"{ROOT}/tmp/{today}-legacy-vs-mp2020-corrections.png")


# %% Python crashes with segfault on correcting the energy of wbm-step-1-24459 due to
# https://github.com/spglib/spglib/issues/194 when using spglib v2.0.{0,1}
cse = df_wbm.computed_structure_entry["wbm-step-1-24459"]
cse = ComputedStructureEntry.from_dict(cse)
mp_compat.process_entry(cse)


# %%
with gzip.open(f"{ROOT}/data/mp/2022-09-18-ppd-mp.pkl.gz", "rb") as zip_file:
    ppd_mp: PatchedPhaseDiagram = pickle.load(zip_file)


# %% calculate e_above_hull for each material
# this loop needs the warnings filter above to not crash Jupyter kernel with logs
# takes ~20 min at 200 it/s for 250k entries in WBM
e_above_hull_key = "e_above_hull_uncorrected_ppd_mp"
assert e_above_hull_key not in df_summary

for entry in tqdm(df_wbm.cse):
    assert entry.entry_id.startswith("wbm-step-")

    e_per_atom = entry.uncorrected_energy_per_atom
    e_hull_per_atom = ppd_mp.get_hull_energy_per_atom(entry.composition)
    e_above_hull = e_per_atom - e_hull_per_atom

    df_summary.at[entry.entry_id, e_above_hull_key] = e_above_hull


# add old + new MP energy corrections to above hull energies
for corrections in ("mp2020", "legacy"):
    df_summary[e_above_hull_key.replace("un", f"{corrections}_")] = (
        df_summary[e_above_hull_key]
        + df_summary[f"e_correction_per_atom_{corrections}"]
    )


# %% calculate formation energies from CSEs wrt MP elemental reference energies
# first make sure source and target dfs have matching indices
assert sum(df_wbm.index != df_summary.index) == 0

e_form_key = "e_form_per_atom_uncorrected_mp_refs"
assert e_form_key not in df_summary

for row in tqdm(df_wbm.itertuples(), total=len(df_wbm)):
    mat_id, cse, formula = row.Index, row.cse, row.formula_from_cse
    assert mat_id == cse.entry_id, f"{mat_id=} != {cse.entry_id=}"
    assert mat_id in df_summary.index, f"{mat_id=} not in df_summary"

    entry_like = dict(composition=formula, energy=cse.uncorrected_energy)
    e_form = get_e_form_per_atom(entry_like)
    e_form_ppd = ppd_mp.get_form_energy_per_atom(cse)

    # make sure the PPD and functional method of calculating formation energy agree
    assert abs(e_form - e_form_ppd) < 1e-7, f"{e_form=} != {e_form_ppd=}"
    df_summary.at[cse.entry_id, e_form_key] = e_form

assert len(df_summary) == sum(
    step_lens
), f"rows were added: {len(df_summary)=} {sum(step_lens)=}"


# add old + new MP energy corrections to formation energies
for corrections in ("mp2020", "legacy"):
    df_summary[e_form_key.replace("un", f"{corrections}_")] = (
        df_summary[e_form_key] + df_summary[f"e_correction_per_atom_{corrections}"]
    )


# %%
df_summary.round(6).to_csv(f"{module_dir}/{today}-wbm-summary.csv")

df_summary = pd.read_csv(f"{module_dir}/2022-10-19-wbm-summary.csv").set_index(
    "material_id"
)


# %% read WBM dataset from disk
df_wbm = pd.read_json(
    f"{module_dir}/2022-10-19-wbm-cses+init-structs.json.bz2"
).set_index("material_id")

df_init_struct = pd.read_json(
    f"{module_dir}/2022-10-19-wbm-init-structs.json.bz2"
).set_index("material_id")

df_wbm["cse"] = [
    ComputedStructureEntry.from_dict(x) for x in tqdm(df_wbm.computed_structure_entry)
]


# %%
df_wbm["init_struct"] = [
    Structure.from_dict(x) if x else None for x in tqdm(df_wbm.initial_structure)
]

wyckoff_col = "wyckoff_spglib"
for idx, struct in tqdm(df_wbm.init_struct.items(), total=len(df_wbm)):
    if struct is None:
        continue
    if not df_wbm.at[idx, wyckoff_col]:
        df_wbm.at[idx, wyckoff_col] = get_aflow_label_from_spglib(struct)


# %% make sure material IDs within each step are consecutive
for step in range(1, 6):
    df = df_summary[df_summary.index.str.startswith(f"wbm-step-{step}-")]
    step_len = step_lens[step - 1]
    assert len(df) == step_len, f"{step=} has {len(df)=}, expected {step_len=}"

    step_counts = list(df.index.str.split("-").str[-1].astype(int))
    assert step_counts == list(
        range(1, step_len + 1)
    ), f"{step=} counts not consecutive"
