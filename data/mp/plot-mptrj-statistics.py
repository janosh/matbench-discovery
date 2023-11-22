import glob

import numpy as np
import pandas as pd
from ase import Atom
from ase.io import read
from matplotlib import pyplot as plt
from tqdm.auto import tqdm

files = glob.glob("./mptrj-gga-ggapu/*.extxyz")

# Get a list of unique elements from all files
unique_elements = set()

# Comment out the following block to speed up analysis
# for file in tqdm(files):
#     atom = read(file, index="0")
#     unique_elements.update(atom.get_chemical_symbols())

for element in [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    89,
    90,
    91,
    92,
    93,
    94,
]:
    unique_elements.update([Atom(element).symbol])

print(unique_elements)

# Create an empty DataFrame with columns for each element
element_columns = list(unique_elements)
columns = (
    [
        "formula",
        "energy_per_atom",
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "stress_yz",
        "stress_xz",
        "stress_xy",
    ]
    + [f"{element}_force" for element in element_columns]
    + [f"{element}_magmom" for element in element_columns]
)

data_all = []

for file in tqdm(files):
    traj = read(file, index=":")
    for atoms in traj:
        stress = atoms.get_stress()
        data_row = {
            "formula": atoms.get_chemical_formula(),
            "energy_per_atom": atoms.get_potential_energy()
            / atoms.get_global_number_of_atoms(),
            "stress_xx": stress[0],
            "stress_yy": stress[1],
            "stress_zz": stress[2],
            "stress_yz": stress[3],
            "stress_xz": stress[4],
            "stress_xy": stress[5],
        }

        # Initialize the element columns with NaN
        data_row.update(
            {f"{element}_force": float("nan") for element in element_columns}
        )
        data_row.update(
            {f"{element}_magmom": float("nan") for element in element_columns}
        )

        # Populate the element columns with force and magmom
        for i, element in enumerate(atoms.get_chemical_symbols()):
            data_row[f"{element}_force"] = np.linalg.norm(atoms.get_forces()[i])
            try:
                data_row[f"{element}_magmom"] = atoms.get_magnetic_moments()[i]
            except Exception:
                continue

        data_all.append(data_row)

df = pd.DataFrame(data_all, columns=columns)

df.to_csv("mptrj-gga-ggapu.csv", index=False)

# Plots


with plt.style.context("default"):
    fig, axes = plt.subplot_mosaic(
        """
        ab
        cd
        """,
        figsize=(6, 3),
        layout="constrained",
    )

    i = "a"
    axes[i].hist(df["energy_per_atom"], bins=100)
    axes[i].set(xlabel="energy [eV/atom]", ylabel="count", yscale="log")

    i = "b"
    force_columns = df.filter(like="_force", axis=1)
    force_series = force_columns.stack()

    axes[i].hist(force_series, bins=100)
    axes[i].set(xlabel=r"force [eV/$\AA$]", ylabel="count", yscale="log")

    i = "c"
    stress_hydrostatic = (df["stress_xx"] + df["stress_yy"] + df["stress_zz"]) / 3.0

    axes[i].hist(stress_hydrostatic, bins=100)
    axes[i].set(
        xlabel="$\\frac{1}{3}\\mathrm{Tr}(\\sigma)$ [eV/$\\AA^3$]",
        ylabel="count",
        yscale="log",
    )

    i = "d"
    magmom_columns = df.filter(like="_magmom", axis=1)
    magmom_series = magmom_columns.stack()

    axes[i].hist(magmom_series, bins=100)
    axes[i].set(xlabel=r"magmom [$\mu_B$]", ylabel="count", yscale="log")

    plt.savefig("mptrj-gga-ggapu-statistics.png", transparent=True)
    plt.show()
