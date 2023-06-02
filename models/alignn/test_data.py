# %% Imports
import pandas as pd
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.data import DATA_FILES, df_wbm


# %% Definitions
DEBUG = False
task_type = "IS2RE"
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "initial_structure"
id_col = "material_id"


# %% Get data
data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
    "IS2RE-debug": f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json-1k-samples.bz2",
}[task_type + ("-debug" if DEBUG else "")]
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

df = pd.read_json(data_path).set_index(id_col)

df[target_col] = df_wbm[target_col]
if task_type == "RS2RE":
    df[input_col] = [x["structure"] for x in df.computed_structure_entry]
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [Structure.from_dict(x) for x in tqdm(df[input_col], disable=None)]


# %% write to ALIGNN format
df[target_col].to_csv("train-data/targets.csv")

for mat_id, struct in tqdm(df[input_col].items(), desc="Saving structures"):
    struct.to(f"{mat_id}.poscar")
