#%% Imports
### ----------------------------------------------------------------------------
import pandas as pd
import os

from pymatgen.io.ase import AseAtomsAdaptor as pase
from pymatgen.core import Structure
from jarvis.core.atoms import Atoms
from sklearn.model_selection import train_test_split
from tqdm import tqdm, trange

from matbench_discovery.data import DATA_FILES
from matbench_discovery.structure import perturb_structure

#%% Definitions
### ----------------------------------------------------------------------------

target_col = "formation_energy_per_atom"
input_col  = "structure"
id_col     = "material_id"
perturb    = 0

#%% Get data
### ----------------------------------------------------------------------------

# Read structures
df_structures = pd.read_json(DATA_FILES.mp_computed_structure_entries).set_index(id_col)
df_structures = pd.Series([Structure.from_dict(s['structure']) for s in tqdm(df_structures['entry'], disable=None)], index=df_structures.index)

# Read energies
df = pd.read_csv(DATA_FILES.mp_energies).set_index(id_col)
df[input_col] = df_structures[df.index]

assert target_col in df

#%% Augment with perturbed structures
### ----------------------------------------------------------------------------

df_aug = df.copy()
structs = df_aug.pop(input_col)
for idx in trange(perturb, desc="Generating perturbed structures"):
    df_aug[input_col] = [perturb_structure(x) for x in structs]
    df = pd.concat([df, df_aug.set_index(f"{x}-aug={idx+1}" for x in df_aug.index)])

del df_aug

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

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=42)

save_data(X_train, y_train, 'data_train')
save_data(X_test , y_test , 'data_test' )

# %%
