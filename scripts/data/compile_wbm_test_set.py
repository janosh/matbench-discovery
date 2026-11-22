"""Download, process and clean the DFT dataset published in Wang et al. 2021.

https://nature.com/articles/s41524-020-00481-6
"""

# %%
import os
from glob import glob

import numpy as np
import pandas as pd
from matscikit import Composition, PhaseDiagram, group_structures
from matscikit.interop import structure_from_pymatgen_dict
from pandas.core.util.hashing import hash_pandas_object
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import WBM_DIR, today
from matbench_discovery.energy import formation_energies
from matbench_discovery.enums import DataFiles, MbdKey
from matbench_discovery.remote.fetch import download_file
from matbench_discovery.wbm import (
    WBM_ELEMENTAL_REF_ENERGIES,
    mp2020_corrected,
    unique_wbm_prototypes,
)


def alphabetical_formula(composition: Composition) -> str:
    """Formula with elements in alphabetical order, e.g. "Ag4 Bi4 O12"."""
    return " ".join(sorted(composition.formula.split()))


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
os.makedirs(f"{WBM_DIR}/raw", exist_ok=True)

drive_download = (
    "https://drive.usercontent.google.com/download?export=download&confirm=t"
)
for step, file_id in google_drive_ids.items():
    file_path = f"{WBM_DIR}/raw/wbm-structures-step-{step}.json.bz2"

    if os.path.isfile(file_path):
        print(f"{file_path} already exists, skipping")
        continue

    print(f"{step=}")
    download_file(file_path, f"{drive_download}&id={file_id}")


# %%
summary_path = f"{WBM_DIR}/raw/wbm-summary.txt"

if not os.path.isfile(summary_path):
    summary_id_file = "1639IFUG7poaDE2uB6aISUOi65ooBwCIg"
    download_file(summary_path, f"{drive_download}&id={summary_id_file}")


# %%
json_paths = sorted(glob(f"{WBM_DIR}/raw/wbm-structures-step-*.json.bz2"))
step_lens = (61848, 52800, 79205, 40328, 23308)
# step 3 has 79,211 initial structures but only 79,205 ComputedStructureEntries
# i.e. 6 extra structures which have missing energy, volume, etc. in the summary file
bad_struct_ids = (70802, 70803, 70825, 70826, 70828, 70829)
# step 5 has 2 missing initial structures: 23166, 23294


assert len(json_paths) == len(step_lens), "Mismatch in WBM steps and JSON files"
wbm_structs_index_checksums = (
    10630821823676988257,
    18360475612623866193,
    10739373004389012550,
    14867548025423706528,
    18198704957443186264,
)

dfs_wbm_structs = locals().get("dfs_wbm_structs", {})

for json_path in json_paths:
    step = int(json_path.split(".json.bz2")[0][-1])
    assert step in range(1, 6)

    if step in dfs_wbm_structs:
        print(f"{step=} {json_path=} already loaded.")
        continue

    print(f"{step=}")
    df_wbm_step = pd.read_json(json_path).T

    # we hash index only for speed
    # could use joblib.hash(df) to hash whole df but it's slow
    checksum = hash_pandas_object(df_wbm_step.index).sum()
    expected = wbm_structs_index_checksums[step - 1]
    assert checksum == expected, (
        f"bad df.index checksum for {step=}, {expected=}, got {checksum=}\n"
        f"\n{json_path=}"
    )

    if step == 3:
        df_wbm_step = df_wbm_step.drop(
            index=[f"step_3_{wbm_id}" for wbm_id in bad_struct_ids]
        )
        # re-index after dropping bad structures to get same indices as summary file
        # where IDs are consecutive, i.e. step_3_70801 is followed by step_3_70802,
        # not step_3_70804, etc.
        # df.index = [f"step_3_{idx + 1}" for idx in range(len(df))]

    step_len = step_lens[step - 1]
    assert len(df_wbm_step) == step_len, (
        f"bad len for {step=}: {len(df_wbm_step)} != {step_len}"
    )
    dfs_wbm_structs[step] = df_wbm_step


# NOTE step 5 is missing 2 initial structures, see nan_init_structs_ids below
assert dict(dfs_wbm_structs[5].isna().sum()) == {"opt": 0, "org": 2}
assert list(dfs_wbm_structs[5].query("org.isna()").index) == [
    "step_5_23165",
    "step_5_23293",
]


# %%
df_wbm = pd.concat(dfs_wbm_structs.values())
del dfs_wbm_structs, df_wbm_step

assert len(df_wbm) == sum(step_lens)


def increment_wbm_material_id(wbm_id: str) -> str:
    """Map step_1_0, step_1_1, ... to wbm-1-1, wbm-1-2, ..."""
    try:
        prefix, step_num, material_num = wbm_id.split("_")
    except ValueError, AttributeError:
        print(f"bad {wbm_id=}")
        return wbm_id

    msg = f"bad {wbm_id=}, {prefix=} {step_num=} {material_num=}"
    assert prefix == "step", msg
    assert step_num.isdigit(), msg
    assert material_num.isdigit(), msg

    return f"wbm-{step_num}-{int(material_num) + 1}"


df_wbm.index = df_wbm.index.map(increment_wbm_material_id)
df_wbm.index.name = Key.mat_id
assert df_wbm.index[0] == "wbm-1-1"
assert df_wbm.index[-1] == "wbm-5-23308"

df_wbm[Key.init_struct] = df_wbm.pop("org")
df_wbm[Key.final_struct] = df_wbm.pop("opt")
assert list(df_wbm.columns) == [Key.init_struct, Key.final_struct]


# %% download WBM ComputedStructureEntries from
# https://archive.materialscloud.org/record/2021.68
mat_cloud_url = "https://archive.materialscloud.org/api/records/qp8qv-gry23/files"

for filename in (
    "README.txt",
    # same file as summary_path above except this doesn't have ID column
    "summary.txt.bz2",
    *(f"step_{step}.json.bz2" for step in range(1, 6)),
):
    file_path = f"{WBM_DIR}/raw/wbm-cse-{filename.lower().replace('_', '-')}"
    if os.path.isfile(file_path):
        print(f"{file_path} already exists, skipping")
        continue

    url = f"{mat_cloud_url}/{filename}/content"
    download_file(file_path, url)


# %%
cse_step_paths = sorted(glob(f"{WBM_DIR}/raw/wbm-cse-step-*.json.bz2"))
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

    df_wbm_step = pd.read_json(json_path)

    step_len = step_lens[step - 1]
    dfs_wbm_cses[step] = df_wbm_step
    assert len(df_wbm_step) == step_len, f"{step=}: {len(df_wbm_step)} != {step_len}"


# %%
df_wbm[Key.computed_structure_entry] = np.concatenate(
    [*dfs_wbm_cses.values()]
).squeeze()
del dfs_wbm_cses, df_wbm_step

for mat_id, cse in df_wbm[Key.computed_structure_entry].items():
    # needed for the MP2020 corrections to classify the entries
    cse["parameters"]["run_type"] = (
        "GGA+U" if cse["parameters"]["is_hubbard"] else "GGA"
    )
    cse["entry_id"] = mat_id
    assert cse["entry_id"].startswith("wbm-")

assert pd.Series(
    cse["parameters"]["run_type"] for cse in tqdm(df_wbm[Key.computed_structure_entry])
).value_counts().to_dict() == {"GGA": 248481, "GGA+U": 9008}

# make sure only 2 materials have missing initial structures with expected IDs
expected_ids_with_no_init_structs = ["wbm-5-23166", "wbm-5-23294"]
actual_ids_with_no_structs = list(df_wbm.query(f"{Key.init_struct}.isna()").index)
assert actual_ids_with_no_structs == expected_ids_with_no_init_structs, (
    f"{actual_ids_with_no_structs=}\n{expected_ids_with_no_init_structs=}"
)
# drop the two materials with missing initial structures
df_wbm = df_wbm.drop(index=expected_ids_with_no_init_structs)


# %% get composition from CSEs
df_wbm["composition_from_cse"] = [
    structure_from_pymatgen_dict(cse["structure"]).composition
    for cse in tqdm(df_wbm[Key.computed_structure_entry])
]

df_wbm["composition_from_final_struct"] = [
    structure_from_pymatgen_dict(struct).composition
    for struct in tqdm(df_wbm[Key.final_struct])
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
for mat_id, row in tqdm(
    df_wbm.sample(n_samples, random_state=0).iterrows(), total=n_samples
):
    struct_final = structure_from_pymatgen_dict(row[Key.final_struct])
    struct_from_cse = structure_from_pymatgen_dict(
        row[Key.computed_structure_entry]["structure"]
    )
    # StructureMatcher's default tolerances and fit criterion (every site within stol);
    # on 1000 random initial/relaxed WBM pairs this agrees with pymatgen's
    # StructureMatcher().fit on 992 (8 more matches: it also searches equivalent
    # reduced bases, so it can match pairs pymatgen rejects)
    assert len(group_structures([struct_final, struct_from_cse])) == 1, (
        f"structure mismatch for {mat_id=}"
    )

    # and check initial and final compositions match
    struct_init = structure_from_pymatgen_dict(row[Key.init_struct])
    assert struct_init.composition == struct_final.composition, (
        f"composition mismatch for {mat_id=}"
    )


del df_wbm[Key.final_struct]


# %% extract alphabetical formula from CSEs (will be used as ground-truth formulas since
# more informative than reduced formulas found in df_summary)
df_wbm["formula_from_cse"] = [
    alphabetical_formula(comp) for comp in df_wbm.pop("composition_from_cse")
]


# %%
col_map = {
    "# comp": Key.formula,
    "nsites": Key.n_sites,
    "vol": "volume",
    "e": MbdKey.dft_energy,
    "e_form": MbdKey.e_form_wbm,
    "e_hull": MbdKey.each_wbm,
    "gap": Key.bandgap_pbe,
    "id": Key.mat_id,
}
# WBM summary was shared twice, once on google drive, once on materials cloud
# download both and check for consistency
df_summary = pd.read_csv(
    f"{WBM_DIR}/raw/wbm-summary.txt", sep="\t", names=list(col_map.values())
).set_index(Key.mat_id)

df_summary_bz2 = pd.read_csv(f"{WBM_DIR}/raw/wbm-cse-summary.txt.bz2", sep="\t").rename(
    columns=col_map
)

# duplicate Ga3Ru2U3 step_3_28147 (1st one is wbm-2-18689) has 0 volume in
# df_summary_bz2 vs 155.41 in df_summary
query_str = "volume > 0 & formula != 'Ga3Ru2U3'"
pd.testing.assert_frame_equal(
    df_summary.reset_index(drop=True).query(query_str),
    df_summary_bz2.reset_index(drop=True).query(query_str),
)

no_id_mask = df_summary.index.isna()
assert sum(no_id_mask) == 6, f"{sum(no_id_mask)=}"
# the 'None' materials have 0 volume, energy, n_sites, bandgap, etc.
assert df_summary[no_id_mask].drop(columns=[Key.formula]).eq(0).to_numpy().all()
assert len(df_summary.query("volume > 0")) == len(df_wbm) + len(
    expected_ids_with_no_init_structs
)
# make sure dropping materials with 0 volume removes exactly 6 materials, the same ones
# listed in bad_struct_ids above
assert all(
    df_summary.reset_index().query("volume == 0").index.to_numpy() - sum(step_lens[:2])
    == bad_struct_ids
)

df_summary.index = df_summary.index.map(increment_wbm_material_id)  # format IDs
# drop materials with id=NaN and missing initial structures
df_summary = df_summary.drop(index=[*expected_ids_with_no_init_structs, float("NaN")])

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
updated_ids: list[str] = []
for mat_id, cse in df_wbm[Key.computed_structure_entry].items():
    entry_id = cse["entry_id"]
    if mat_id != entry_id:
        print(f"{mat_id=} != {entry_id=}, updating entry_id to mat_id")
        cse["entry_id"] = mat_id
        updated_ids += [entry_id]

print(f"total mismatching CSE IDs fixed: {len(updated_ids):,}")

# Normalize source formulas before comparing reduced and unreduced compositions.
df_summary[Key.formula] = [
    alphabetical_formula(Composition(formula)) for formula in df_summary[Key.formula]
]


# %% write initial structures and computed structure entries to compressed json
for fname, cols in (
    ("computed-structure-entries", [Key.computed_structure_entry]),
    ("init-structs", [Key.init_struct]),
    (
        "computed-structure-entries+init-structs",
        [Key.init_struct, Key.computed_structure_entry],
    ),
):
    all_cols = ["formula_from_cse", *cols]
    df_wbm[all_cols].reset_index().to_json(f"{WBM_DIR}/{today}-wbm-{fname}.json.bz2")


del df_wbm[Key.init_struct]


# %% df_summary and df_wbm formulas differ because summary formulas are reduced while
# df_wbm formulas are not (e.g. Ac6 U2 vs Ac3 U1 in summary). unreduced is more
# informative so we use it.
assert sum(df_summary[Key.formula] != df_wbm.formula_from_cse) == 114_273
assert sum(df_summary[Key.formula] == df_wbm.formula_from_cse) == 143_214

df_summary[Key.formula] = df_wbm.formula_from_cse


# fix bad energy which is 0 in df_summary but a more realistic -63.68 in CSE
df_summary.loc["wbm-2-18689", MbdKey.dft_energy] = df_wbm.loc["wbm-2-18689"][
    Key.computed_structure_entry
]["energy"]

# NOTE in ComputedEntry dicts, "energy" is the uncorrected energy (pymatgen's
# entry.energy property is uncorrected + correction)


# %% compare summary energies with CSE energies
df_summary[f"{MbdKey.dft_energy}_from_cse"] = [
    cse["energy"] for cse in tqdm(df_wbm[Key.computed_structure_entry])
]

# check CSE and summary energies are consistent, only exceeding 0.1 eV difference twice
diff_e_cse_e_summary = (
    df_summary[MbdKey.dft_energy] - df_summary.uncorrected_energy_from_cse
)
assert diff_e_cse_e_summary.max() < 0.15
assert sum(diff_e_cse_e_summary > 0.1) == 2


# %% remove suspicious formation energy outliers
e_form_cutoff = 5
n_too_stable = sum(df_summary[MbdKey.e_form_wbm] < -e_form_cutoff)
print(f"{n_too_stable = }")  # n_too_stable = 502
n_too_unstable = sum(df_summary[MbdKey.e_form_wbm] > e_form_cutoff)
print(f"{n_too_unstable = }")  # n_too_unstable = 22


assert len(df_summary) == len(df_wbm) == 257_487

query_str = f"{-e_form_cutoff} < {MbdKey.e_form_wbm} < {e_form_cutoff}"
dropped_ids = sorted(set(df_summary.index) - set(df_summary.query(query_str).index))
assert len(dropped_ids) == 502 + 22
assert dropped_ids[:3] == ["wbm-1-12142", "wbm-1-12143", "wbm-1-12144"]
assert dropped_ids[-3:] == ["wbm-5-9121", "wbm-5-9211", "wbm-5-934"]

df_summary = df_summary.query(query_str)
df_wbm = df_wbm.loc[df_summary.index]


# make sure we dropped the expected number 524 of materials
assert len(df_summary) == len(df_wbm) == 256_963  # 257_487 - 502 - 22


# %%
for mat_id, cse in df_wbm[Key.computed_structure_entry].items():
    assert mat_id == cse["entry_id"], f"{mat_id} != {cse['entry_id']}"

wbm_cses = list(df_wbm[Key.computed_structure_entry])
# raw WBM ComputedStructureEntries have no energy corrections applied:
assert all(cse["correction"] == 0 for cse in wbm_cses)
# summary and CSE n_sites match
assert all(df_summary.n_sites == [len(cse["structure"]["sites"]) for cse in wbm_cses])

# current MP2020 corrections (pymatgen's defaults: strict_anions="require_bound" since
# pymatgen#3803, an anion correction needs oxidation state <= -1), the scheme model
# predictions are corrected with. The legacy pre-2020 MaterialsProjectCompatibility
# column e_correction_per_atom_mp_legacy is carried over from the published summary
# below (no current implementation), so older models trained on the legacy scheme
# (e.g. MEGNet) can still swap corrections.
df_wbm_corr = mp2020_corrected(wbm_cses)
pd.testing.assert_index_equal(df_wbm_corr.index, df_summary.index, check_names=False)

n_corrected = (df_wbm_corr.correction_per_atom != 0).sum()
assert n_corrected == 98_711, f"{n_corrected=:,}, expected 98,711"

df_summary["e_correction_per_atom_mp2020"] = df_wbm_corr.correction_per_atom.to_numpy()

avg_energy_corr = df_summary.e_correction_per_atom_mp2020.mean().round(4)
assert avg_energy_corr == -0.1062, f"{avg_energy_corr=:.4}, expected -0.1062"


# %% calculate e_above_hull for each material
if MbdKey.each_true in df_summary:
    raise KeyError(f"{MbdKey.each_true!s} already in {df_summary.columns=}")

if mismatched_ids := [
    mat_id
    for mat_id, cse in df_wbm[Key.computed_structure_entry].items()
    if mat_id != cse["entry_id"]
]:
    raise ValueError(f"entry_id != material_id for {mismatched_ids[:10]=}")
# the hull spans every WBM composition: a ValueError here means a WBM entry has an
# element no MP entry has
wbm_compositions = list(df_wbm_corr.composition)
corrected_energies = df_wbm_corr.energy_per_atom.to_numpy()
raw_energies = df_wbm_corr.uncorrected_energy.to_numpy()
del df_wbm, wbm_cses, df_wbm_corr


# %% convex hull of all MP entries (MP2020-corrected energies), solved per query as a
# linear program by matscikit: ~1 s to build and query all ~257k WBM entries vs ~6 min
# for the pymatgen PatchedPhaseDiagram loop it replaces (plus ~5 min to build the PPD)
# MP entries are recorrected with the same scheme as the WBM entries: their stored
# corrections came from the MP database build and differ for 329 of 154,718 entries
df_mp_entries = pd.read_json(DataFiles.mp_computed_structure_entries.path)
df_mp_corr = mp2020_corrected(list(df_mp_entries.entry))
assert len(df_mp_corr) == len(df_mp_entries), (
    f"{len(df_mp_corr)=} {len(df_mp_entries)=}"
)
mp_hull = PhaseDiagram(
    list(df_mp_corr.composition), df_mp_corr.energy_per_atom.to_numpy()
)
del df_mp_entries, df_mp_corr


df_summary[MbdKey.each_true] = mp_hull.e_above_hull(
    wbm_compositions, corrected_energies
)


# %% calculate formation energies from CSEs wrt MP elemental reference energies
e_forms_raw = formation_energies(
    wbm_compositions,
    raw_energies,
    WBM_ELEMENTAL_REF_ENERGIES,
)
# make sure the phase diagram and standalone reference-energy methods of calculating
# formation energy agree
e_forms_hull = (
    np.asarray(mp_hull.formation_energy(wbm_compositions, corrected_energies))
    - df_summary["e_correction_per_atom_mp2020"].to_numpy()
)
e_form_diffs = pd.Series(e_forms_raw - e_forms_hull, index=df_summary.index)
if (bad_e_forms := e_form_diffs[e_form_diffs.abs() >= 1e-4]).any():
    raise ValueError(f"formation energies disagree by >= 1e-4 eV/atom: {bad_e_forms}")
df_summary[MbdKey.e_form_raw] = e_forms_raw


df_summary[MbdKey.e_form_raw.replace("uncorrected", "mp2020_corrected")] = (
    df_summary[MbdKey.e_form_raw] + df_summary["e_correction_per_atom_mp2020"]
)


# %% Keep the structural descriptors that define the published benchmark. Its spglib
# labels can differ from current assignments even after converting Pearson symbols.
archived_cols = [
    "e_correction_per_atom_mp_legacy",
    MbdKey.init_protostructure_spglib,
    MbdKey.protostructure_spglib,
    "site_stats_fingerprint_init_final_norm_diff",
]
archived = pd.read_csv(
    DataFiles.wbm_summary_2023.path,
    index_col=Key.mat_id,
    usecols=[Key.mat_id, *archived_cols, MbdKey.uniq_proto],
)
pd.testing.assert_index_equal(df_summary.index, archived.index)
df_summary = df_summary.join(archived[archived_cols], validate="one_to_one")
n_legacy = (df_summary.e_correction_per_atom_mp_legacy != 0).sum()
assert n_legacy == 39_591, f"{n_legacy=:,} expected 39,591"

# Rebuild cohort membership from relaxed labels and original WBM hull distances.
df_mp = pd.read_csv(DataFiles.mp_energies.path, usecols=[MbdKey.protostructure_spglib])
df_summary[MbdKey.uniq_proto] = unique_wbm_prototypes(
    df_summary, df_mp[MbdKey.protostructure_spglib]
)
pd.testing.assert_series_equal(
    df_summary[MbdKey.uniq_proto], archived[MbdKey.uniq_proto]
)


# %% write final summary data to disk (yeah!)
schema = pd.read_csv(
    DataFiles.wbm_summary_2023.path, nrows=0, index_col=Key.mat_id
).columns
pd.testing.assert_index_equal(df_summary.columns.sort_values(), schema.sort_values())
df_summary[schema].round(6).to_csv(f"{WBM_DIR}/{today}-wbm-summary.csv.gz")
