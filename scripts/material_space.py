# %% Imports
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# pip install umap-learn
import umap
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery.data import DATA_FILES


# %% Definitions
basename = "material_space"
data_dir = "../data"

filename_features_mp = os.path.join(
    data_dir, "mp/2023-02-07-mp-computed-structure-entries-features.csv.bz2"
)
filename_features_wbm = os.path.join(
    data_dir, "wbm/2022-10-19-wbm-computed-structure-entries-features.csv.bz2"
)
filename_umap = os.path.join(data_dir, f"{basename}.npz")
filename_plot = f"{basename}.png"

n_neighbors = 750


# %% Material Structure and Composition Featurizer
# "A critical examination of robustness and generalizability of machine learning prediction
#  of materials properties"](https://www.nature.com/articles/s41524-023-01012-9) by
# Kangming Li, Brian DeCost, Kamal Choudhary, Michael Greenwood, and Jason Hattrick-Simpers.


def to_unitcell(structure):
    [site.to_unit_cell(in_place=True) for site in structure.sites]
    return structure


def featurize_dataframe(df_in, col_id="structure", ignore_errors=True, chunksize=30):
    """Featurize a dataframe using Matminter Structure featurizer

    Parameters
    ----------
    df : Pandas.DataFrame
        DataFrame with a column named "structure"

    Returns:
    -------
    A DataFrame containing 273 features (columns)

    """
    # For featurization
    from matminer.featurizers.base import MultipleFeaturizer
    from matminer.featurizers.composition import (
        ElementProperty,
        IonProperty,
        Stoichiometry,
        ValenceOrbital,
    )
    from matminer.featurizers.structure import (
        ChemicalOrdering,
        MaximumPackingEfficiency,
        SiteStatsFingerprint,
        StructuralHeterogeneity,
        StructureComposition,
    )

    if isinstance(df_in, pd.Series):
        df = df_in.to_frame()
    else:
        df = df_in
    df[col_id] = df[col_id].apply(to_unitcell)

    # 128 structural feature
    struc_feat = [
        SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
        SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
        StructuralHeterogeneity(),
        MaximumPackingEfficiency(),
        ChemicalOrdering(),
    ]
    # 145 compositional features
    compo_feat = [
        StructureComposition(Stoichiometry()),
        StructureComposition(ElementProperty.from_preset("magpie")),
        StructureComposition(ValenceOrbital(props=["frac"])),
        StructureComposition(IonProperty(fast=True)),
    ]
    featurizer = MultipleFeaturizer(struc_feat + compo_feat)
    # Set the chunksize used for Pool.map parallelisation
    featurizer.set_chunksize(chunksize=chunksize)
    featurizer.fit(df[col_id])
    X = featurizer.featurize_dataframe(
        df=df, col_id=col_id, ignore_errors=ignore_errors
    )
    # check failed entries
    print("Featurization completed.")
    failed = np.any(pd.isnull(X.iloc[:, df.shape[1] :]), axis=1)
    if np.sum(failed) > 0:
        print(f"Number failed: {np.sum(failed)}/{len(failed)}")
    return X


# %%
def featurize_file(
    filename,
    computed_structure=False,
    input_col="initial_structure",
    id_col="material_id",
):
    df_in = pd.read_json(filename).set_index(id_col)

    if computed_structure:
        df_in[input_col] = [
            Structure.from_dict(x["structure"])
            for x in tqdm(
                df_in[input_col], leave=False, desc="Converting to PyMatgen Structure"
            )
        ]
    else:
        df_in[input_col] = [
            Structure.from_dict(x)
            for x in tqdm(
                df_in[input_col], leave=False, desc="Converting to PyMatgen Structure"
            )
        ]

    df_features = featurize_dataframe(df_in[input_col], col_id=input_col)

    return df_features.drop(input_col, axis=1)


# %%
def select_features(df_features, threshold=0.95):
    corr_matrix = df_features.corr(method="pearson", numeric_only=False).abs()

    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

    # Find features with correlation greater than our threshold
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]

    # Drop features
    return df_features.drop(to_drop, axis=1)


def import_features(filename, id_col="material_id"):
    df = pd.read_csv(filename)
    df = df.set_index(id_col)
    return df


# %%
def plot_scatter(u, y, filename=None):
    plt.scatter(u[:, 0], u[:, 1], c=y, cmap="Spectral", s=5)
    plt.colorbar(boundaries=np.arange(min(y), max(y) + 2) - 0.5).set_ticks(
        np.arange(min(y), max(y) + 1)
    )
    if filename is not None:
        plt.savefig(filename)
    plt.show()


# %% Compute features and export to CSV
if not os.path.exists(filename_features_mp):
    print("Computing matminer mp features...")

    df_features = featurize_file(
        DATA_FILES.mp_computed_structure_entries,
        input_col="entry",
        computed_structure=True,
    )
    df_features.to_csv(filename_features_mp)

if not os.path.exists(filename_features_wbm):
    print("Computing matminer wbm features...")

    df_features = featurize_file(DATA_FILES.wbm_initial_structures)
    df_features.to_csv(filename_features_wbm)


# %%
if not os.path.exists(filename_umap):
    print("Computing UMAP...")

    df_in = import_features(filename_features_wbm)
    df_mp = import_features(filename_features_mp)
    df = pd.concat((df_in, df_mp))

    # Drop all rows containing NaN values
    df_in = df_in.dropna(axis=0)
    df_mp = df_mp.dropna(axis=0)

    # Combined data frame
    df = pd.concat((df_in, df_mp))

    # Drop highly correlated features
    df_mp = select_features(df_mp, threshold=0.95)
    df_in = df_in[df_mp.columns]
    df = df[df_mp.columns]

    # %% Create labels
    y = [re.sub(r"^wbm-(\d+)-\d+$", r"\1", id) for id in df.index]
    y = [re.sub(r"^(?:mp|mvc)-\d+$", r"0", id) for id in y]
    y = pd.Series(y, index=df.index).astype(int).to_numpy()

    # Train only on MP data
    reducer = umap.UMAP(random_state=42, low_memory=False, n_neighbors=n_neighbors)
    reducer.fit(df_mp)

    # Transform everything
    u = reducer.transform(df)

    np.savez(filename_umap, u=u, y=y)


# %%
if not os.path.exists(filename_plot):
    print("Plotting material space...")

    u = np.load(filename_umap)["u"]
    y = np.load(filename_umap)["y"]

    plot_scatter(u, y, filename=filename_plot)
