"""Test EquFlash model on matbench-discovery IS2RE task."""

import argparse
import json
import multiprocessing as mp
import os
import pickle
from typing import Any

import pandas as pd
import tqdm
from GGNN.common.calculator import UCalculator
from GGNN.preprocessing.atoms_to_graphs import AtomsToGraphs
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from relaxation.ml_relaxation import ml_relax
from relaxation.optimizable import batch_to_atoms
from torch_geometric.data import Batch

from matbench_discovery.data import as_dict_handler
from matbench_discovery.hpc import df_slurm_chunk


def load_pickle(file_path: str) -> Any:  # noqa: ANN401
    """Load a single pickle file."""
    with open(file_path, "rb") as file:
        return pickle.load(file)  # noqa: S301 safe internal data


def load_multiple_pickles(file_paths: list[str]) -> list[Any]:
    """Load multiple pickle files in parallel."""
    n_proc = min(32, mp.cpu_count())
    with mp.Pool(processes=n_proc) as pool:
        return pool.map(load_pickle, file_paths)


def split_df_by_max_atoms(
    df: pd.DataFrame,
    max_atoms: int,
    input_col: str = "initial_structure",
) -> list[pd.DataFrame]:
    """Split dataframe so each batch has at most max_atoms atoms."""
    batches = []
    current_rows = []
    current_atoms = 0

    for idx, row in df.iterrows():
        struct = Structure.from_dict(row[input_col])
        n_atoms = len(struct)

        if n_atoms > max_atoms:
            raise ValueError(
                f"{idx} has {n_atoms} atoms, larger than max_atoms={max_atoms}"
            )

        if current_rows and current_atoms + n_atoms > max_atoms:
            batches.append(df.loc[current_rows])
            current_rows = []
            current_atoms = 0

        current_rows.append(idx)
        current_atoms += n_atoms

    if current_rows:
        batches.append(df.loc[current_rows])

    return batches


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, help="Model checkpoint path")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--rank", type=int, default=0)
    parser.add_argument("--fmax", type=float, default=0.02)
    parser.add_argument("--worldsize", type=int, default=1)
    parser.add_argument("--max-atoms", type=int, default=1000)
    parser.add_argument(
        "--init-structs-dir", required=True, help="WBM initial structures directory"
    )
    parser.add_argument(
        "--wbm-metadata-file",
        required=True,
        help="WBM structures metadata file",
    )
    return parser.parse_args()


if __name__ == "__main__":
    """Run IS2RE evaluation on matbench-discovery WBM test set."""
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    metadata = {
        "checkpoint": os.path.abspath(args.checkpoint),
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
    df = df_slurm_chunk(df, args.worldsize, args.rank + 1)  # rank is 0-based

    input_col = "initial_structure"
    batch_lists = split_df_by_max_atoms(df, args.max_atoms, input_col=input_col)

    calc = UCalculator(checkpoint_path=args.checkpoint, cpu=False)
    trainer = calc.trainer

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

        for mat_id, atoms, _ntraj in zip(
            mat_idx, atoms_list, n_traj.numpy(), strict=True
        ):
            unwrapped = getattr(atoms, "atoms", atoms)
            relaxed_struct = AseAtomsAdaptor.get_structure(unwrapped)
            relax_results[mat_id] = {
                "energy": atoms.get_total_energy(),
                "structure": relaxed_struct,
            }

    df_out = pd.DataFrame(relax_results).T.add_prefix("mlff_")
    df_out.index.name = Key.mat_id
    df_metadatas = pd.read_csv(args.wbm_metadata_file)
    df_metadatas = df_metadatas.set_index("material_id").reindex(df_out.index)
    if df_metadatas.isna().all(axis=1).any():
        missing = df_metadatas[df_metadatas.isna().all(axis=1)].index.tolist()
        raise ValueError(f"Missing metadata for material IDs: {missing[:5]}...")
    df_merge = pd.concat([df_metadatas, df_out], axis=1)
    df_merge["ggnn_e_per_form"] = (
        df_merge["mlff_energy"]
        + df_merge["correction"]
        - df_merge["formation_ref_energy"]
    ) / df_merge["num_atoms"]
    structure_outfile = (
        f"{args.out}/{args.rank:03d}_{args.worldsize}_structures.jsonl.gz"
    )

    df_struct = df_merge["mlff_structure"].reset_index()

    df_struct.to_json(
        structure_outfile,
        orient="records",
        lines=True,
        compression="gzip",
        default_handler=as_dict_handler,
    )

    df_merge.drop(columns=["mlff_structure"]).reset_index().to_json(
        output_jsonfile_name,
        default_handler=as_dict_handler,
        compression="gzip",
    )
