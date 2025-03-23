"""
The original relaxation trajectories shared by Microsoft Research contained xyz=(0,0,0)
for all atom positions.

This script written on 2025-03-23 converts the pickled relaxation directories that MSR
shared upon second request into standard Matbench Discovery JSON Lines format for
geometry optimization analysis.
Original submission date: 2024-12-19

The original trajectories can be found at:
https://figshare.com/s/a629acbf3bed6a04b3ce?file=53060504
"""

import pickle

import pandas as pd
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import ROOT, today

for model_key, model_name in {
    "mattersim_1M": "mattersim-v1-1M",
    "mattersim_5M": "mattersim-v1-5M",
}.items():
    pickle_path = f"{ROOT}/models/mattersim/{model_name}/{model_key}.json"
    with open(pickle_path, mode="rb") as file:
        data = pickle.load(file)  # noqa: S301

    struct_col = f"{model_key}_structure"
    df_geo_opt = pd.DataFrame(
        {**traj[-1].info, struct_col: AseAtomsAdaptor.get_structure(traj[-1])}
        for traj in tqdm(data.values())
    )

    df_geo_opt = df_geo_opt.set_index("material_id")
    df_geo_opt = df_geo_opt.drop(columns=["structure_index"])

    json_path = f"{ROOT}/models/mattersim/{model_name}/{today}-wbm-geo-opt.jsonl.gz"
    df_geo_opt.reset_index().to_json(
        json_path,
        orient="records",
        lines=True,
        default_handler=lambda x: x.as_dict() if hasattr(x, "as_dict") else x,
    )

    # load first 5 back from JSON for sanity check
    df_test = pd.read_json(json_path, lines=True, nrows=5)
    print(f"{df_test=}")
