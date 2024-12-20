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

# Path to the tar.gz archive
module_dir = os.path.dirname(__file__)
tgz_path = f"{module_dir}/2024-12-19-wbm-relax-traj-mattersim-v1-5M.tgz"

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
pmg_json_path = "mattersim-v1-5M-relaxed-structures.json.gz"
pd.Series(structures_dict).to_json(pmg_json_path)


# %% reload from disk
srs_relaxed = pd.read_json(pmg_json_path, typ="series")
print(f"{len(srs_relaxed)=:,}")
