"""Test EquFlash model on matbench-discovery IS2RE task."""
# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
# "torch==2.8.0+cu126",
# "torch-geometric==2.6.1",
# "numpy==1.26.0",
# "scikit-learn==1.7.2",
# "spglib==2.6.0",
# "e3nn==0.5.6",
# "ase==3.26.0",
# "pymatgen==2025.10.7",
# "pymatviz>=0.16.0",
# "flashTP_e3nn==0.1.0",
# "fairchem-core==1.10.0",
# "lmdb==1.6.2",
# "submitit==1.5.3",
# "matbench-discovery==1.3.1",
# ]
#
# [tool.uv.sources]
# flashTP_e3nn = { git = "https://github.com/SNU-ARC/flashTP" }
# matbench-discovery = { path = "../../", editable = true }
# ///

import argparse
import json
import multiprocessing as mp
import os
import pickle
from typing import Any

import numpy as np
import pandas as pd
import torch
import tqdm
from fairchem.core.common.utils import update_config
from GGNN.preprocessing import AtomsToGraphs
from GGNN.trainer.utrainer import OCPTrainer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import SiteCollection, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from relaxation.ml_relaxation import ml_relax
from relaxation.optimizable import batch_to_atoms
from torch_geometric.data import Batch


def split_df(df: pd.DataFrame, batch_size: int) -> list[pd.DataFrame]:
    """Split dataframe into batches."""
    return [df[idx : idx + batch_size] for idx in range(0, len(df), batch_size)]


def compute_rmsd(
    struct_pred: SiteCollection, struct_og: SiteCollection
) -> tuple[float | None, float | None]:
    """Compute RMSD between predicted and ground truth structures."""
    structure_matcher = StructureMatcher(stol=1.0, scale=False)
    result = structure_matcher.get_rms_dist(struct_pred, struct_og)
    return result if result else (None, None)


def load_trainer_from_ckpt(checkpoint_path: str) -> OCPTrainer:
    """Load EquFlash trainer from checkpoint."""
    from fairchem.core.common.utils import setup_imports

    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    config = checkpoint["config"]

    if "model_attributes" in config:
        config["model_attributes"]["name"] = config.pop("model")
        config["model"] = config["model_attributes"]

    config["model"]["otf_graph"] = True
    config = update_config(config)
    config["checkpoint"] = str(checkpoint_path)
    del config["dataset"]["src"]

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
    """Serialize MSONable objects to dict. Non-serializable objects become None."""
    try:
        return obj.as_dict()
    except AttributeError:
        return None


def load_pickle(file_path: str) -> Any:
    """Load a single pickle file."""
    with open(file_path, "rb") as file:
        return pickle.load(file)  # noqa: S301 safe internal data


def load_multiple_pickles(file_paths: list[str]) -> list[Any]:
    """Load multiple pickle files in parallel."""
    n_proc = min(32, mp.cpu_count())
    with mp.Pool(processes=n_proc) as pool:
        return pool.map(load_pickle, file_paths)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--ckpt", required=True, help="Model checkpoint path")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--rank", type=int, default=0)
    parser.add_argument("--fmax", type=float, default=0.02)
    parser.add_argument("--worldsize", type=int, default=1)
    parser.add_argument("--batchsize", type=int, default=32)
    parser.add_argument(
        "--init-structs-dir", required=True, help="WBM initial structures directory"
    )
    parser.add_argument(
        "--relaxed-structs-file",
        required=True,
        help="WBM relaxed structures pickle file",
    )
    return parser.parse_args()


if __name__ == "__main__":
    """Run IS2RE evaluation on matbench-discovery WBM test set."""
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    with open(args.relaxed_structs_file, "rb") as file:
        structs_wbm = pickle.load(file)  # noqa: S301 safe internal data

    metadata = {
        "ckpt": os.path.abspath(args.ckpt),
        "worldsize": args.worldsize,
        "fmax": args.fmax,
    }
    output_jsonfile_name = f"{args.out}/{args.rank:03d}_{args.worldsize}.json.gz"
    metadata_json_path = f"{args.out}/metadata_{args.rank:03d}_{args.worldsize}.json"

    if os.path.isfile(metadata_json_path):
        with open(metadata_json_path) as file:
            metadata_old = json.load(file)
        for key, val in metadata.items():
            if val != metadata_old[key]:
                raise FileExistsError(
                    "Output folder already exists with conflicting metadata. "
                    "Please change --out."
                )
    else:
        with open(metadata_json_path, "w") as file:
            json.dump(metadata, file)

    if os.path.isfile(output_jsonfile_name):
        print("Results already exist!")
        raise SystemExit(0)

    file_lists = [
        f"{args.init_structs_dir}/{fname}"
        for fname in sorted(os.listdir(args.init_structs_dir))
    ]
    dfs = load_multiple_pickles(file_lists)
    df = pd.concat(dfs, axis=0)
    df = np.array_split(df, args.worldsize)[args.rank]

    input_col = "initial_structure"
    batch_lists = split_df(df, args.batchsize)

    trainer = load_trainer_from_ckpt(checkpoint_path=args.ckpt)

    a2g = AtomsToGraphs(r_edges=False)
    relax_results = {}
    for batch_df in tqdm.tqdm(batch_lists):
        atoms_list = [
            a2g.convert(Structure.from_dict(atoms_dict).to_ase_atoms())
            for atoms_dict in batch_df[input_col].to_numpy().tolist()
        ]
        mat_idx = batch_df.index.tolist()
        for mat_id, struct in zip(mat_idx, atoms_list, strict=True):
            struct.sid = mat_id
        batch = Batch.from_data_list(atoms_list)
        n_traj, batch = ml_relax(batch, trainer, 500, args.fmax, relax_cell=True)
        atoms_list = batch_to_atoms(
            batch,
            {"forces": batch.forces, "energy": batch.energy, "stress": batch.stress},
        )

        for mat_id, atoms, ntraj in zip(
            mat_idx, atoms_list, n_traj.numpy(), strict=True
        ):
            unwrapped = getattr(atoms, "atoms", atoms)
            relaxed_struct = AseAtomsAdaptor.get_structure(unwrapped)
            rmsd, max_dist = compute_rmsd(relaxed_struct, structs_wbm[mat_id])
            relax_results[mat_id] = {
                "energy": atoms.get_total_energy(),
                "rmsd": rmsd,
                "max_dist": max_dist,
                "n_traj": ntraj,
            }

    df_out = pd.DataFrame(relax_results).T.add_prefix("mlff_")
    df_out.index.name = Key.mat_id

    df_out.reset_index().to_json(output_jsonfile_name, default_handler=as_dict_handler)
