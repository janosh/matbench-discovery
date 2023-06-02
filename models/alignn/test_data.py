#%% Imports
### ----------------------------------------------------------------------------
import pandas as pd
import os

from pymatgen.io.ase import AseAtomsAdaptor as pase
from pymatgen.core import Structure
from jarvis.core.atoms import Atoms
from tqdm import tqdm

from matbench_discovery import DEBUG, ROOT
from matbench_discovery.data import DATA_FILES, df_wbm

#%% Definitions
### ----------------------------------------------------------------------------

DEBUG      = False
task_type  = "IS2RE"
target_col = "e_form_per_atom_mp2020_corrected"
input_col  = "initial_structure"
id_col     = "material_id"
perturb    = 0

#%% Get data
### ----------------------------------------------------------------------------

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

#%% Convert ase to jarvis
### ----------------------------------------------------------------------------

def ase_to_atoms(ase_atoms="", cartesian=True):
    """Convert ase structure to Atoms."""
    return Atoms(
        lattice_mat = ase_atoms.get_cell(),
        elements    = ase_atoms.get_chemical_symbols(),
        coords      = ase_atoms.get_positions(),
        cartesian   = cartesian,
    )

#%% Generic function to export data in ALIGNN compatible format
### ----------------------------------------------------------------------------

def save_data(X, y, basename):

    if not os.path.exists(basename):
        os.mkdir(basename)

    with open(os.path.join(basename, 'id_prop.csv'), 'w') as f:
        for i, x in enumerate(X):
            filename_poscar = f'POSCAR-{i}.vasp'
            # Convert structure
            x_ase  = pase.get_atoms(x)
            x_atom = ase_to_atoms(x_ase)
            # Export structure to poscar file
            x_atom.write_poscar(os.path.join(basename, filename_poscar))
            # Write filename and target value to CSV file
            f.write('%s,%6f\n' % (filename_poscar, y[i]))

#%% Export data
### ----------------------------------------------------------------------------

X = list(df[input_col])
y = list(df[target_col])

save_data(X, y, 'data_test_wbm' )

# %%
