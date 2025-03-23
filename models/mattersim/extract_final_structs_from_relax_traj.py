import os
import tarfile
from io import BytesIO
from typing import Any

import pandas as pd

# unused import needed to parse string repr of Atoms with eval() below
from ase.atoms import Atoms  # noqa: F401
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import Model

# Path to the tar.gz archive
module_dir = os.path.dirname(__file__)
model_key = Model.mattersim_v1_5m.key
tgz_path = f"{module_dir}/{model_key}/2024-12-19-wbm-relax-traj.tgz"

# Dictionary to store structures
structures_dict: dict[str, dict[str, Any]] = {}
ase_adaptor = AseAtomsAdaptor()

# Open and process the tar.gz archive
with tarfile.open(tgz_path, "r:gz") as tar:
    # Get list of CSV files in the archive
    extracted_dir = tar.getnames()[0]

    # Get all files in the archive
    all_members = tar.getmembers()
    # Filter for CSV files
    csv_members = [m for m in all_members if m.name.endswith(".csv.gz")]

    # Process each trajectory file
    for member in tqdm(csv_members, desc="Processing trajectory files"):
        # Extract material_id from filename
        file = os.path.basename(member.name)
        idx = int(file.split("_")[-1].split(".")[0])
        mat_id = df_wbm.iloc[idx][Key.mat_id]

        # Read the compressed CSV file directly from tar archive
        csv_gz_file = tar.extractfile(member)
        if csv_gz_file is None:
            raise ValueError(f"No CSV file found in {member.name}")
        df_i = pd.read_csv(BytesIO(csv_gz_file.read()))

        # Get the final structure (last row) from the 'Atoms' column
        final_atoms = eval(df_i["atoms"].iloc[-1])  # noqa: S307

        # Convert ASE Atoms to Pymatgen Structure
        structure = ase_adaptor.get_structure(final_atoms)

        # Store in dictionary with material_id as key
        structures_dict[mat_id] = structure.as_dict()


# Save to compressed JSON
pmg_json_path = f"{module_dir}/{model_key}/2024-12-19-wbm-geo-opt.jsonl.gz"
df_structs = pd.Series(structures_dict).to_frame()
df_structs.index.name = Key.mat_id
df_structs.columns = [f"{model_key}_{Key.structure}"]
df_structs.reset_index().to_json(pmg_json_path)
print(f"Saved {len(df_structs):,} structures to {pmg_json_path}")


# %% reload from disk
df_from_disk = pd.read_json(pmg_json_path, typ="series").to_frame()
print(f"{len(df_from_disk)=:,}")
