"""This script computes the UMAP projections of matminer features for the WBM and MP
datasets. The resulting projections are saved to a .npz file.

Needs pip install umap-learn.

Uses Material Structure and Composition Featurizer
"A critical examination of robustness and generalizability of machine learning
prediction of materials properties"
https://www.nature.com/articles/s41524-023-01012-9
"""

# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymatviz as pmv
import umap
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR, PDF_FIGS, WBM_DIR
from matbench_discovery.data import DataFiles

__author__ = "Philipp Benner, Janosh Riebesell"
__date__ = "2023-11-28"


# %%
mp_matminer_feat_path = f"{MP_DIR}/mp-final-structures-matminer-features.csv.bz2"
wbm_matminer_feat_path = f"{WBM_DIR}/wbm-initial-structures-matminer-features.csv.bz2"
# default n_umap_neighbors=15 results resemble PCA, i.e. a single blob of points.
# The larger n_neighbors, the more MP point clusters move apart. Those islands may
# correspond to some material classes (but we didn't test this).
n_umap_neighbors = 750


# %%
def featurize_dataframe(
    df_in: pd.DataFrame | pd.Series,
    *,
    struct_col: str = "structure",
    ignore_errors: bool = True,
    chunk_size: int = 30,
) -> pd.DataFrame:
    """Featurize a dataframe of structures using matminer Featurizers.

    Args:
        df_in (pd.DataFrame): DataFrame with a column named "structure"
        struct_col (str, optional): Name of column containing structures.
            Defaults to "structure".
        ignore_errors (bool, optional): Whether to ignore errors. Defaults to True.
        chunk_size (int, optional): Chunk size for parallelization. Defaults to 30.

    Returns:
        pd.DataFrame: DataFrame with 273 columns containing matminer features
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

    df_in = df_in.to_frame() if isinstance(df_in, pd.Series) else df_in
    for struct in df_in[struct_col]:
        for site in struct:
            site.to_unit_cell(in_place=True)

    # 128 structural feature
    struct_feat = [
        SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
        SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
        StructuralHeterogeneity(),
        MaximumPackingEfficiency(),
        ChemicalOrdering(),
    ]
    # 145 compositional features
    comp_feat = [
        StructureComposition(Stoichiometry()),
        StructureComposition(ElementProperty.from_preset("magpie")),
        StructureComposition(ValenceOrbital(props=["frac"])),
        StructureComposition(IonProperty(fast=True)),
    ]
    featurizer = MultipleFeaturizer([*struct_feat, *comp_feat])
    # Set the chunk size used for Pool.map parallelization
    featurizer.set_chunksize(chunksize=chunk_size)
    featurizer.fit(df_in[struct_col])
    matminer_features = featurizer.featurize_dataframe(
        df=df_in, col_id=struct_col, ignore_errors=ignore_errors
    )
    print("Featurization complete")
    # check failed entries
    failed = np.any(pd.isna(matminer_features.iloc[:, df_in.shape[1] :]), axis=1)
    if sum(failed) > 0:
        print(f"Number failed: {sum(failed)} / {len(failed)}")
    return matminer_features


# %%
def featurize_file(filename: str, struct_col: str = Key.init_struct) -> pd.DataFrame:
    """Featurize pymatgen Structures in a file with matminer."""
    df_in = pd.read_json(filename).set_index(Key.mat_id)

    # ComputedStructureEntry dicts have a "structure" key, if that's missing it's a
    # Structure dict
    df_in[struct_col] = [
        Structure.from_dict(x.get("structure", x))
        for x in tqdm(df_in[struct_col], leave=False, desc="Hydrate structures")
    ]

    df_features = featurize_dataframe(df_in[struct_col], struct_col=struct_col)

    return df_features.drop(columns=struct_col)


# %%
def features_to_drop(df_in: pd.DataFrame, threshold: float = 0.95) -> list[str]:
    """Get column names of features with correlation greater than threshold to drop."""
    corr_matrix = df_in.corr(method="pearson", numeric_only=False).abs()

    # select upper triangle of correlation matrix
    upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

    # return features with correlation greater than the threshold
    return [col for col in upper_tri if any(upper_tri[col] > threshold)]


# %% Compute matminer features for MP and WBM, then export to CSV
if not os.path.isfile(mp_matminer_feat_path):
    df_mp_feats = featurize_file(
        DataFiles.mp_computed_structure_entries.path, struct_col="entry"
    )
    df_mp_feats.to_csv(mp_matminer_feat_path)

if not os.path.isfile(wbm_matminer_feat_path):
    df_wbm_feats = featurize_file(DataFiles.wbm_initial_structures.path)
    df_wbm_feats.to_csv(wbm_matminer_feat_path)


# %% Compute UMAP projection of matminer features
umap_out_path = f"{WBM_DIR}/umap/2d-umap-projections.csv.bz2"
if not os.path.isfile(umap_out_path):
    df_mp = pd.read_csv(mp_matminer_feat_path).set_index(Key.mat_id)
    df_wbm = pd.read_csv(wbm_matminer_feat_path).set_index(Key.mat_id)

    # Drop all rows containing NaN values
    df_wbm = df_wbm.dropna(axis="index")
    df_mp = df_mp.dropna(axis="index")

    # Drop highly correlated features
    cols_to_drop = features_to_drop(df_mp, threshold=0.95)
    df_mp = df_mp.drop(columns=cols_to_drop)
    df_wbm = df_wbm.drop(columns=cols_to_drop)

    # Combined data frame
    df_all = pd.concat((df_wbm, df_mp))
    df_all.index = df_all.index.str.replace(r"^(mp|mvc)-", "wbm-0")
    df_all["wbm_step"] = df_all.index.str.split("-").str[1].astype(int)

    # train only on MP data
    reducer = umap.UMAP(
        n_neighbors=n_umap_neighbors, n_components=2, random_state=42, low_memory=False
    )
    reducer.fit(df_mp)

    # transform everything
    umap_points = reducer.transform(df_all)
    umap_cols = [f"UMAP {idx + 1}" for idx in range(umap_points.shape[1])]
    df_umap = pd.DataFrame(umap_points, index=df_all.index, columns=umap_cols)

    df_umap.to_csv(umap_out_path)


# %% Plot UMAP 2d projection
df_umap = pd.read_csv(umap_out_path).set_index("wbm_step")

umap_cols = list(df_umap)
if umap_cols != ["UMAP 1", "UMAP 2"]:
    raise ValueError(f"Unexpected {umap_cols=}")
min_step, max_step = df_umap.index.min(), df_umap.index.max()
ax = df_umap.plot.scatter(
    *umap_cols, c=df_umap.index, cmap="Spectral", s=5, figsize=(6, 4), colorbar=False
)
cbar = ax.figure.colorbar(
    ax.collections[0],
    boundaries=np.arange(min_step, max_step + 2) - 0.5,
    ticks=range(min_step, max_step + 1),
)
cbar.ax.set_title("WBM step (0 = MP)", rotation=90, y=0.5, x=3, va="center")


# %%
plt.tight_layout()
pmv.save_fig(ax, f"{PDF_FIGS}/wbm-final-struct-matminer-features-2d-umap.png", dpi=300)
