# %% Imports
import os

import pandas as pd
from pymatgen.core import Structure
from sklearn.model_selection import train_test_split
from tqdm import tqdm, trange

from matbench_discovery.data import DATA_FILES
from matbench_discovery.structure import perturb_structure

__author__ = "Philipp Benner"
__date__ = "2023-06-02"


# %%
target_col = "formation_energy_per_atom"
input_col = "structure"
id_col = "material_id"
n_perturb = 0


# %% load structures
df_cse = pd.read_json(DATA_FILES.mp_computed_structure_entries).set_index(id_col)
df_cse[input_col] = [
    Structure.from_dict(cse[input_col]) for cse in tqdm(df_cse.entry, disable=None)
]

# load energies
df = pd.read_csv(DATA_FILES.mp_energies).set_index(id_col)
df[input_col] = df_cse[input_col]
assert target_col in df


# %% augment with randomly perturbed structures
df_aug = df.copy()
structs = df_aug.pop(input_col)
for idx in trange(n_perturb, desc="Generating perturbed structures"):
    df_aug[input_col] = [perturb_structure(x) for x in structs]
    df = pd.concat([df, df_aug.set_index(f"{x}-aug={idx+1}" for x in df_aug.index)])

del df_aug


# %% export data
X_train, X_test, y_train, y_test = train_test_split(
    df[input_col], df[target_col], test_size=0.05, random_state=42
)

for samples, targets, label in ((X_train, y_train, "train"), (X_test, y_test, "test")):
    out_dir = f"{label}-data"
    os.makedirs(out_dir, exist_ok=True)

    targets.to_csv(f"{out_dir}/targets.csv")

    struct: Structure
    for mat_id, struct in tqdm(
        samples.items(), desc="Saving structures", total=len(samples)
    ):
        struct.to(f"{out_dir}/{mat_id}.poscar", fmt="POSCAR")
