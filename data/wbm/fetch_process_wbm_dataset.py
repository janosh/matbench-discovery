# %%
import gzip
import os
import pickle
import urllib.request
import warnings
from glob import glob

import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from pymatgen.core import Composition, Structure
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility,
    MaterialsProjectCompatibility,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz import density_scatter
from pymatviz.utils import save_fig
from tqdm import tqdm

from matbench_discovery import SITE_FIGS, today
from matbench_discovery.data import DATA_FILES
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.plots import pio

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
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

assert pio.templates.default == "plotly_dark+global"


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
os.makedirs(f"{module_dir}/raw", exist_ok=True)

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
# step 3 has 79,211 initial structures but only 79,205 ComputedStructureEntries
# i.e. 6 extra structures which have missing energy, volume, etc. in the summary file
bad_struct_ids = (70802, 70803, 70825, 70826, 70828, 70829)
# step 5 has 2 missing initial structures: 23166, 23294


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
        # df.index = [f"step_3_{idx + 1}" for idx in range(len(df))]

    step_len = step_lens[step - 1]
    assert len(df) == step_len, f"bad len for {step=}: {len(df)} != {step_len}"
    dfs_wbm_structs[step] = df


# NOTE step 5 is missing 2 initial structures, see nan_init_structs_ids below
assert dict(dfs_wbm_structs[5].isna().sum()) == {"opt": 0, "org": 2}
assert list(dfs_wbm_structs[5].query("org.isna()").index) == [
    "step_5_23165",
    "step_5_23293",
]


# %%
df_wbm = pd.concat(dfs_wbm_structs.values())

assert len(df_wbm) == sum(step_lens)


def increment_wbm_material_id(wbm_id: str) -> str:
    """Map step_1_0, step_1_1, ... to wbm-1-1, wbm-1-2, ..."""
    try:
        prefix, step_num, material_num = wbm_id.split("_")
    except (ValueError, AttributeError):
        print(f"bad {wbm_id=}")
        return wbm_id

    msg = f"bad {wbm_id=}, {prefix=} {step_num=} {material_num=}"
    assert prefix == "step", msg
    assert step_num.isdigit(), msg
    assert material_num.isdigit(), msg

    return f"wbm-{step_num}-{int(material_num) + 1}"


df_wbm.index = df_wbm.index.map(increment_wbm_material_id)
df_wbm.index.name = "material_id"
assert df_wbm.index[0] == "wbm-1-1"
assert df_wbm.index[-1] == "wbm-5-23308"

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

    try:
        url = f"{mat_cloud_url}&filename={filename}"
        urllib.request.urlretrieve(url, file_path)
    except urllib.error.HTTPError as exc:
        print(f"failed to download {url=}: {exc}")
        continue


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
df_wbm["computed_structure_entry"] = np.concatenate([*dfs_wbm_cses.values()]).squeeze()

for mat_id, cse in df_wbm.computed_structure_entry.items():
    # needed to ensure MaterialsProjectCompatibility can process the entries
    cse["parameters"]["run_type"] = (
        "GGA+U" if cse["parameters"]["is_hubbard"] else "GGA"
    )
    cse["entry_id"] = mat_id
    assert cse["entry_id"].startswith("wbm-")

assert pd.Series(
    cse["parameters"]["run_type"] for cse in tqdm(df_wbm.computed_structure_entry)
).value_counts().to_dict() == {"GGA": 248481, "GGA+U": 9008}

# make sure only 2 materials have missing initial structures with expected IDs
nan_init_structs_ids = ["wbm-5-23166", "wbm-5-23294"]
assert list(df_wbm.query("initial_structure.isna()").index) == nan_init_structs_ids
# drop the two materials with missing initial structures
df_wbm = df_wbm.drop(index=nan_init_structs_ids)


# %% get composition from CSEs
df_wbm["composition_from_cse"] = [
    ComputedStructureEntry.from_dict(cse).composition
    for cse in tqdm(df_wbm.computed_structure_entry)
]

df_wbm["composition_from_final_struct"] = [
    Structure.from_dict(struct).composition for struct in tqdm(df_wbm.final_structure)
]

# all but 1 composition matches between CSE and final structure
# mismatching ID: wbm-1-37977 which becomes equal on reduction:
# CSE Comp: Ag4 Bi4 O12
# final structure Comp: Ag16 Bi16 O48
df_mismatch = df_wbm.query("composition_from_cse != composition_from_final_struct")
assert len(df_mismatch) == 1
assert df_mismatch.index[0] == "wbm-1-37977"
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


# %% extract alphabetical formula from CSEs (will be used as ground-truth formulas since
# more informative than reduced formulas found in df_summary)
df_wbm["formula_from_cse"] = [
    x.alphabetical_formula for x in df_wbm.pop("composition_from_cse")
]


# %%
col_map = {
    "# comp": "formula",
    "nsites": "n_sites",
    "vol": "volume",
    "e": "uncorrected_energy",
    "e_form": "e_form_per_atom_wbm",
    "e_hull": "e_above_hull_wbm",
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

# duplicate Ga3Ru2U3 step_3_28147 (1st one is wbm-2-18689) has 0 volume in
# df_summary_bz2 vs 155.41 in df_summary
query_str = "volume > 0 & formula != 'Ga3Ru2U3'"
pd.testing.assert_frame_equal(
    df_summary.reset_index(drop=True).query(query_str),
    df_summary_bz2.reset_index(drop=True).query(query_str),
)


assert sum(no_id_mask := df_summary.index.isna()) == 6, f"{sum(no_id_mask)=}"
# the 'None' materials have 0 volume, energy, n_sites, bandgap, etc.
assert all(df_summary[no_id_mask].drop(columns=["formula"]) == 0)
assert len(df_summary.query("volume > 0")) == len(df_wbm) + len(nan_init_structs_ids)
# make sure dropping materials with 0 volume removes exactly 6 materials, the same ones
# listed in bad_struct_ids above
assert all(
    df_summary.reset_index().query("volume == 0").index.to_numpy() - sum(step_lens[:2])
    == bad_struct_ids
)

df_summary.index = df_summary.index.map(increment_wbm_material_id)  # format IDs
# drop materials with id=NaN and missing initial structures
df_summary = df_summary.drop(index=[*nan_init_structs_ids, float("NaN")])

# the 8403 material IDs in step 3 with final number larger than any of the ones in
# bad_struct_ids are now misaligned between df_summary and df_wbm
# the IDs in df_summary are consecutive while the IDs in df_wbm skip over the numbers in
# bad_struct_ids. we fix this with fix_bad_struct_index_mismatch() by mapping the IDs in
# df_wbm to the ones in df_summary so that both indices become consecutive.
assert sum(df_summary.index != df_wbm.index) == 8403
assert {*df_summary.index} - {*df_wbm.index} == {
    "wbm-3-70803",
    "wbm-3-70804",
    "wbm-3-70826",
    "wbm-3-70827",
    "wbm-3-70829",
    "wbm-3-70830",
}


def fix_bad_struct_index_mismatch(material_id: str) -> str:
    """Decrement material IDs in step 3 by the number of IDs with smaller final number
    in bad_struct_ids. This should fix the index mismatch between df_summary and df_wbm.
    """
    _, step_num, mat_num = material_id.split("-")
    step_num, mat_num = int(step_num), int(mat_num)

    if step_num == 3:
        mat_num -= sum(mat_num > idx + 1 for idx in bad_struct_ids)

    return f"wbm-{step_num}-{mat_num}"


# don't accidentally apply the fix twice
if sum(df_summary.index != df_wbm.index) != 0:
    df_wbm.index = df_wbm.index.map(fix_bad_struct_index_mismatch)

# check that the index mismatch is fixed
assert sum(df_summary.index != df_wbm.index) == 0

# update ComputedStructureEntry entry_ids to match material_ids
for mat_id, cse in df_wbm.computed_structure_entry.items():
    entry_id = cse["entry_id"]
    if mat_id != entry_id:
        print(f"{mat_id=} != {entry_id=}, updating entry_id to mat_id")
        cse["entry_id"] = mat_id


# sort formulas alphabetically
df_summary["alph_formula"] = [
    Composition(x).alphabetical_formula for x in df_summary.formula
]
# alphabetical formula and original formula differ due to spaces, number 1 after element
# symbols (FeO vs Fe1 O1), and element order (FeO vs OFe)
assert sum(df_summary.alph_formula != df_summary.formula) == 257_483

df_summary["formula"] = df_summary.pop("alph_formula")


# %% write initial structures and computed structure entries to compressed json
for fname, cols in (
    ("computed-structure-entries", ["computed_structure_entry"]),
    ("init-structs", ["initial_structure"]),
    (
        "computed-structure-entries+init-structs",
        ["initial_structure", "computed_structure_entry"],
    ),
):
    cols = ["formula_from_cse", *cols]
    df_wbm[cols].reset_index().to_json(f"{module_dir}/{today}-wbm-{fname}.json.bz2")


# %%
# df_summary and df_wbm formulas differ because summary formulas are reduced while
# df_wbm formulas are not (e.g. Ac6 U2 vs Ac3 U1 in summary). unreduced is more
# informative so we use it.
assert sum(df_summary.formula != df_wbm.formula_from_cse) == 114_273
assert sum(df_summary.formula == df_wbm.formula_from_cse) == 143_214

df_summary.formula = df_wbm.formula_from_cse


# fix bad energy which is 0 in df_summary but a more realistic -63.68 in CSE
df_summary.loc["wbm-2-18689", "uncorrected_energy"] = df_wbm.loc[
    "wbm-2-18689"
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
assert diff_e_cse_e_summary.max() < 0.15
assert sum(diff_e_cse_e_summary > 0.1) == 2

density_scatter(df_summary.uncorrected_energy, df_summary.uncorrected_energy_from_cse)


# %% remove suspicious formation energy outliers
e_form_cutoff = 5
n_too_stable = sum(df_summary.e_form_per_atom_wbm < -e_form_cutoff)
print(f"{n_too_stable = }")  # n_too_stable = 502
n_too_unstable = sum(df_summary.e_form_per_atom_wbm > e_form_cutoff)
print(f"{n_too_unstable = }")  # n_too_unstable = 22

fig = df_summary.hist(
    x="e_form_per_atom_wbm", backend="plotly", log_y=True, range_x=[-5.5, 5.5]
)
fig_compressed = False
fig.add_vline(x=e_form_cutoff, line=dict(dash="dash"))
fig.add_vline(x=-e_form_cutoff, line=dict(dash="dash"))
fig.add_annotation(
    **dict(x=0, y=1, yref="paper", yshift=20),
    text=f"<b>dataset cropped to within +/- {e_form_cutoff} eV/atom</b>",
    showarrow=False,
)
x_axis_title = "WBM uncorrected formation energy (eV/atom)"
fig.update_layout(xaxis_title=x_axis_title, margin=dict(l=10, r=10, t=40, b=10))
fig.update_yaxes(fixedrange=True)  # disable zooming y-axis
fig.show(
    config=dict(
        modeBarButtonsToRemove=["lasso2d", "select2d", "autoScale2d", "toImage"],
        displaylogo=False,
    )
)


# %%
# no need to store all 250k x values in plot, leads to 1.7 MB file, subsample every 10th
# point is enough to see the distribution, round to 3 decimal places to reduce file size
if not fig_compressed:
    fig_compressed = True
    fig.data[0].x = [round(x, 3) for x in fig.data[0].x[::10]]

img_path = f"{SITE_FIGS}/hist-wbm-e-form-per-atom"
save_fig(fig, f"{img_path}.svelte")
# recommended to upload SVG to vecta.io/nano for compression
# save_fig(fig, f"{img_path}.svg", width=800, height=300)


# %%
assert len(df_summary) == len(df_wbm) == 257_487

query_str = f"{-e_form_cutoff} < e_form_per_atom_wbm < {e_form_cutoff}"
dropped_ids = sorted(set(df_summary.index) - set(df_summary.query(query_str).index))
assert len(dropped_ids) == 502 + 22
assert dropped_ids[:3] == "wbm-1-12142 wbm-1-12143 wbm-1-12144".split()
assert dropped_ids[-3:] == "wbm-5-9121 wbm-5-9211 wbm-5-934".split()

df_summary = df_summary.query(query_str)
df_wbm = df_wbm.loc[df_summary.index]


# make sure we dropped the expected number 524 of materials
assert len(df_summary) == len(df_wbm) == 257_487 - 502 - 22


# %%
for mat_id, cse in df_wbm.computed_structure_entry.items():
    assert mat_id == cse["entry_id"], f"{mat_id} != {cse['entry_id']}"

df_wbm["cse"] = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(df_wbm.computed_structure_entry)
]
# raw WBM ComputedStructureEntries have no energy corrections applied:
assert all(cse.uncorrected_energy == cse.energy for cse in df_wbm.cse)
# summary and CSE n_sites match
assert all(df_summary.n_sites == [len(cse.structure) for cse in df_wbm.cse])


# entries are corrected in-place by default so we apply legacy corrections first
# and then leave the new corrections in place below
# having both old and new corrections allows updating predictions from older models
# like MEGNet that were trained on MP release prior to new corrections by subtracting
# old corrections and adding the new ones
entries_old_corr = MaterialsProjectCompatibility().process_entries(
    df_wbm.cse, clean=True, verbose=True
)
assert len(entries_old_corr) == len(df_wbm), f"{len(entries_old_corr)=} {len(df_wbm)=}"

# extract legacy MP energy corrections to df_megnet
df_wbm["e_correction_per_atom_mp_legacy"] = [
    cse.correction_per_atom for cse in df_wbm.cse
]

# clean up legacy corrections and apply new corrections
entries_new_corr = MaterialsProject2020Compatibility().process_entries(
    df_wbm.cse, clean=True, verbose=True
)
assert len(entries_new_corr) == len(df_wbm), f"{len(entries_new_corr)=} {len(df_wbm)=}"

n_corrected = sum(cse.uncorrected_energy != cse.energy for cse in df_wbm.cse)
assert n_corrected == 100_930, f"{n_corrected=} expected 100,930"

df_summary["e_correction_per_atom_mp2020"] = [
    cse.correction_per_atom for cse in df_wbm.cse
]

assert df_summary.e_correction_per_atom_mp2020.mean().round(4) == -0.1069


# %%
with gzip.open(DATA_FILES.mp_patched_phase_diagram, "rb") as zip_file:
    ppd_mp: PatchedPhaseDiagram = pickle.load(zip_file)


# %% calculate e_above_hull for each material
# this loop needs above warnings.filterwarnings() to not crash Jupyter kernel with logs
# takes ~20 min at 200 it/s for 250k entries in WBM
each_col = "e_above_hull_mp2020_corrected_ppd_mp"
assert each_col not in df_summary

for mat_id, cse in tqdm(df_wbm.cse.items(), total=len(df_wbm)):
    assert mat_id == cse.entry_id, f"{mat_id=} != {cse.entry_id=}"
    assert cse.entry_id in df_summary.index, f"{cse.entry_id=} not in df_summary"

    e_above_hull = ppd_mp.get_e_above_hull(cse, allow_negative=True)

    df_summary.loc[cse.entry_id, each_col] = e_above_hull


# %% calculate formation energies from CSEs wrt MP elemental reference energies
# first make sure source and target dfs have matching indices
assert sum(df_wbm.index != df_summary.index) == 0

e_form_col = "e_form_per_atom_uncorrected"

for row in tqdm(df_wbm.itertuples(), total=len(df_wbm)):
    mat_id, cse, formula = row.Index, row.cse, row.formula_from_cse
    assert mat_id == cse.entry_id, f"{mat_id=} != {cse.entry_id=}"
    assert mat_id in df_summary.index, f"{mat_id=} not in df_summary"

    entry_like = dict(composition=formula, energy=cse.uncorrected_energy)
    e_form = get_e_form_per_atom(entry_like)
    e_form_ppd = ppd_mp.get_form_energy_per_atom(cse) - cse.correction_per_atom

    # make sure the PPD.get_e_form_per_atom() and standalone get_e_form_per_atom()
    # method of calculating formation energy agree
    assert (
        abs(e_form - e_form_ppd) < 1e-4
    ), f"{mat_id}: {e_form=:.3} != {e_form_ppd=:.3} (diff={e_form - e_form_ppd:.3}))"
    df_summary.loc[cse.entry_id, e_form_col] = e_form


df_summary[e_form_col.replace("uncorrected", "mp2020_corrected")] = (
    df_summary[e_form_col] + df_summary["e_correction_per_atom_mp2020"]
)


# %%
try:
    from aviary.wren.utils import get_aflow_label_from_spglib

    wyckoff_col = "wyckoff_spglib"
    if wyckoff_col not in df_wbm:
        df_summary[wyckoff_col] = None

    for idx, struct in tqdm(df_wbm.initial_structure.items(), total=len(df_wbm)):
        if not pd.isna(df_summary.loc[idx, wyckoff_col]):
            continue  # Aflow label already computed
        try:
            struct = Structure.from_dict(struct)
            df_summary.loc[idx, wyckoff_col] = get_aflow_label_from_spglib(struct)
        except Exception as exc:
            print(f"{idx=} {exc=}")

    assert df_summary[wyckoff_col].isna().sum() == 0
except ImportError:
    print("aviary not installed, skipping Wyckoff label generation")
except Exception as exception:
    print(f"Generating Aflow labels raised {exception=}")


# %%
fingerprints_path = f"{module_dir}/site-stats.json.gz"
suggest = "not found, run scripts/compute_struct_fingerprints.py to generate"
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"
try:
    df_fp = pd.read_json(fingerprints_path).set_index("material_id")
    df_summary[fp_diff_col] = df_fp[fp_diff_col]
except FileNotFoundError:
    print(f"{fingerprints_path=} {suggest}")
except KeyError:
    print(f"{fingerprints_path=} does not contain {fp_diff_col=}")


# %% write final summary data to disk (yeah!)
df_summary.round(6).to_csv(f"{module_dir}/{today}-wbm-summary.csv")


# %% only here to load data for later inspection
if False:
    df_summary = pd.read_csv(f"{module_dir}/2022-10-19-wbm-summary.csv.gz").set_index(
        "material_id"
    )
    df_wbm = pd.read_json(
        f"{module_dir}/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2"
    ).set_index("material_id")

    df_wbm["cse"] = [
        ComputedStructureEntry.from_dict(dct)
        for dct in tqdm(df_wbm.computed_structure_entry)
    ]
