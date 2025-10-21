import torch
from pathlib import Path
from copy import deepcopy
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
import tqdm
from typing import Any


def split_df(df, batch_size):
    return [df[i : i + batch_size] for i in range(0, len(df), batch_size)]


# Example usage:

import argparse
import os
import pickle
import pandas as pd
from pymatgen.core import Structure
import numpy as np
import json

import multiprocessing as mp
import pickle
import os
from fairchem.core.common.utils import (
    update_config,
)
from relaxation.ml_relaxation import ml_relax
from GGNN.trainer.utrainer import OCPTrainer

BASEDIR = "/home/gpu1/zetta/aixsim/1_umlff/datasets/matbench-discovery/1.0.0"
with open(f"{BASEDIR}/wbm/2025-03-26-wbm-relaxed-structures.pkl", "rb") as f:
    structs_wbm = pickle.load(f)
from pymatgen.analysis.structure_matcher import StructureMatcher
from relaxation.ase_utils import batch_to_atoms


def compute_rmsd(struct_pred, struct_og):
    structure_matcher = StructureMatcher(stol=1.0, scale=False)
    rmsd, max_dist = structure_matcher.get_rms_dist(struct_pred, struct_og) or (
        None,
        None,
    )
    return (rmsd, max_dist)


def load_trainer_from_ckpt(checkpoint_path):
    checkpoint = torch.load(
        checkpoint_path,
        map_location=torch.device("cpu"),
        weights_only=False,
    )
    config = checkpoint["config"]

    if "model_attributes" in config:
        config["model_attributes"]["name"] = config.pop("model")
        config["model"] = config["model_attributes"]

        # Calculate the edge indices on the fly
    config["model"]["otf_graph"] = True
    ### backwards compatability with OCP v<2.0
    config = update_config(config)
    config["checkpoint"] = str(checkpoint_path)
    del config["dataset"]["src"]
    from fairchem.core.common.utils import setup_imports

    setup_imports()
    trainer = OCPTrainer(
        task=config.get("task", {}),
        model=config["model"],
        dataset=[config["dataset"]],
        outputs=config["outputs"],
        loss_functions=config["loss_functions"],
        evaluation_metrics=config["evaluation_metrics"],
        optimizer=config["optim"],
        identifier="",
        slurm=config.get("slurm", {}),
        local_rank=config.get("local_rank", 0),
        is_debug=config.get("is_debug", True),
        cpu=False,
        amp=config.get("amp", False),
        inference_only=True,
    )
    trainer.load_checkpoint(checkpoint_path, checkpoint, inference_only=True)

    return trainer


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    serialize Python classes with as_dict(). Warning: Objects without a as_dict() method
    are replaced with None in the serialized data.
    """
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories


def load_pickle(file_path):
    """Function to load a single pickle file."""
    with open(file_path, "rb") as f:
        return pickle.load(f)


def load_multiple_pickles(file_paths):
    with mp.Pool(processes=32) as pool:
        results = pool.map(load_pickle, file_paths)
    return results


def parse_args():
    parser = argparse.ArgumentParser("batched matbench evaluate")
    parser.add_argument("--ckpt", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--rank", type=int, default=0)
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--worldsize", type=int, default=1)
    parser.add_argument("--batchsize", type=int, default=32)
    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    metadata = {
        "ckpt": os.path.abspath(args.ckpt),
        "worldsize": args.worldsize,
        "fmax": args.fmax,
    }
    output_jsonfile_name = os.path.join(
        args.out, f"{args.rank:03d}_{args.worldsize}.json.gz"
    )
    metadata_json_path = os.path.join(
        args.out, f"metadata_{args.rank:03d}_{args.worldsize}.json"
    )
    if os.path.isfile(metadata_json_path):
        with open(metadata_json_path, "r") as f:
            metadata_old = json.load(f)
        for k, v in metadata.items():
            if v != metadata_old[k]:
                assert (), "output folder already exists. Please change --out"
    else:
        with open(metadata_json_path, "w") as f:
            json.dump(metadata, f)

    if os.path.isfile(output_jsonfile_name):
        print("results already exists!")
        return

    # load list of snapshots
    wbm_path = "/home/gpu1/zetta/aixsim/1_umlff/datasets/2022-10-19-wbm-init-structs/"
    dfs = []

    file_lists = os.listdir(wbm_path)
    file_lists = [os.path.join(wbm_path, fname) for fname in file_lists]
    file_lists.sort()
    dfs = load_multiple_pickles(file_lists)
    df = pd.concat(dfs, axis=0)
    df = np.array_split(df, args.worldsize)[args.rank]

    # split atoms w.r.t worldsize and rank
    input_col = "initial_structure"
    structs = df[input_col].map(Structure.from_dict).to_dict()
    batch_size = args.batchsize

    batch_lists = split_df(df, batch_size)

    trainer = load_trainer_from_ckpt(checkpoint_path=args.ckpt)
    from GGNN.preprocessing.atoms_to_graphs import AtomsToGraphs
    from torch_geometric.data import Batch

    a2g = AtomsToGraphs(r_edges=False)
    relax_results = {}
    for batch_df in tqdm.tqdm(batch_lists):

        atoms_list = [
            a2g.convert(Structure.from_dict(a).to_ase_atoms())
            for a in batch_df[input_col].values.tolist()
        ]
        mat_idx = batch_df.index.tolist()
        for mat_id, stct in zip(mat_idx, atoms_list):
            stct.sid = mat_id
        batch = Batch.from_data_list(atoms_list)
        n_traj, batch = ml_relax(batch, trainer, 500, args.fmax, relax_cell=True)

        atoms_list = batch_to_atoms(
            batch,
            {"forces": batch.forces, "energy": batch.energy, "stress": batch.stress},
        )
        # optimizable.reset(atoms_list)

        for mat_id, atoms, ntraj in zip(mat_idx, atoms_list, n_traj.numpy()):
            relaxed_struct = AseAtomsAdaptor.get_structure(
                getattr(atoms, "atoms", atoms)
            )
            (rmsd, max_dist) = compute_rmsd(relaxed_struct, structs_wbm[mat_id]) or (
                None,
                None,
            )
            relax_results[mat_id] = {
                "energy": atoms.get_total_energy(),
                "rmsd": rmsd,
                "max_dist": max_dist,
                "n_traj": ntraj,
            }

    df_out = pd.DataFrame(relax_results).T.add_prefix("mlff_")
    df_out.index.name = Key.mat_id

    df_out.reset_index().to_json(output_jsonfile_name, default_handler=as_dict_handler)


if __name__ == "__main__":
    main()
