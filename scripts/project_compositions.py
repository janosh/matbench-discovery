"""Compute t-SNE or UMAP projections of WBM and MP compositions."""

# %%
import os
from datetime import UTC, datetime
from typing import Any, Literal

import numpy as np
import pandas as pd
from pymatgen.core import Composition
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import DATA_DIR, timestamp
from matbench_discovery.data import DataFiles
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-28"


data_name = "mp"  # which data to project
projection_type: Literal["tsne", "umap"] = "tsne"  # which projection method to use
out_dim = 2  # number of dimensions to project to
one_hot_dim = 112  # number of elements to use for one-hot encoding
job_name = f"{data_name}-{projection_type}-{out_dim}d"

out_dir = f"{DATA_DIR}/{data_name}/{projection_type}"
os.makedirs(out_dir, exist_ok=True)

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="6:0:0",
)

data_path = {"wbm": DataFiles.wbm_summary.path, "mp": DataFiles.mp_energies.path}[
    data_name
]
print(f"{data_path=}")
print(f"{out_dim=}")
print(f"{projection_type=}")
start_time = datetime.now(tz=UTC)
print(f"job {job_name} started at {timestamp}")
df_in = pd.read_csv(data_path, na_filter=False).set_index(Key.mat_id)


def metric(
    x: np.ndarray,
    y: np.ndarray,
    err_weight: float = 3,
    split_dim: int = one_hot_dim,
) -> float:
    """Custom metric for t-SNE/UMAP that weights the error dimension higher by a factor
    of err_weight than the composition dimensions.
    """
    x_comp, x_err = np.split(x, [split_dim])
    y_comp, y_err = np.split(y, [split_dim])
    return np.linalg.norm(x_comp - y_comp) + err_weight * np.linalg.norm(x_err - y_err)


if projection_type == "tsne":
    from sklearn.manifold import TSNE

    projector = TSNE(
        n_components=out_dim, random_state=0, n_iter=250, n_iter_without_progress=50
    )
    out_cols = [f"{out_dim}d t-SNE {idx + 1}" for idx in range(out_dim)]
elif projection_type == "umap":
    from umap import UMAP

    # TODO this execution path is untested (was never run yet)
    projector = UMAP(n_components=out_dim, random_state=0, metric=metric)
    out_cols = [f"{out_dim}d UMAP {idx + 1}" for idx in range(out_dim)]

identity = np.eye(one_hot_dim)


def sum_one_hot_elem(formula: str) -> np.ndarray[Any, np.int64]:
    """Return sum of one-hot encoded elements in weighted by amount in composition."""
    return sum(identity[el.Z - 1] * amt for el, amt in Composition(formula).items())


one_hot_encoding = np.array(
    [sum_one_hot_elem(formula) for formula in tqdm(df_in[Key.formula])]
)

projections = projector.fit_transform(one_hot_encoding)

df_in[out_cols] = projections

out_path = f"{out_dir}/one-hot-{one_hot_dim}-composition-{out_dim}d.csv.gz"
df_in[out_cols].to_csv(out_path)

print(f"Wrote projections to {out_path!r}")
end_time = datetime.now(tz=UTC)
print(
    f"Job finished at {end_time:%Y-%m-%d %H:%M:%S} and took "
    f"{(end_time - start_time).seconds} sec"
)
