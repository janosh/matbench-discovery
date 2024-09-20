import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import typer
import wandb
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from orb_models.forcefield.calculator import ORBCalculator
from orb_models.forcefield.pretrained import ORB_PRETRAINED_MODELS
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter

"""
pip install git+https://github.com/janosh/matbench-discovery.git@5c8601a

# Required for valid data paths.
git clone https://github.com/janosh/matbench-discovery.git
git checkout 5c8601a
cd matbench-discovery
pip install -e .
"""


torch.set_float32_matmul_precision("high")

app = typer.Typer(pretty_exceptions_enable=False, no_args_is_help=True)

FILTERS = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}

OPTIMIZERS = {
    "FIRE": FIRE,
    "LBFGS": LBFGS,
}

PREDICTED_ENERGY_COL = "orb_energy"


@app.command()
def main(
    model_name: str = "orb-v1",  # Or orb-v1-mptrj-only
    ase_optimizer: str = typer.Option("FIRE", help="ASE optimizer to use"),
    ase_filter: str = typer.Option("frechet", help="ASE filter to use"),
    device: str = typer.Option(
        "cuda" if torch.cuda.is_available() else "cpu", help="Device to use"
    ),
    out_dir: Path = typer.Option("matbench_eval", help="Output directory"),  # noqa: B008
    max_steps: int = typer.Option(500, help="Max steps"),
    force_max: float = typer.Option(0.05, help="Max force"),
    cell_opt: bool = typer.Option(True, help="Optimize cell"),  # noqa: FBT001, FBT003
    limit: int | None = typer.Option(None, help="Debug mode, only use 100 samples"),
    shard: int | None = typer.Option(None, help="Shard the data"),
    total_shards: int | None = typer.Option(None, help="Total number of shards"),
) -> None:
    """Run ORB relaxation on the WBM dataset.

    Produces (possibly sharded) compressed JSON files with relaxed structures
    and energies. These are then aggregated and evaluated in the
    join_predictions script.

    Raises:
        ValueError: If total_shards and shard are not both None or both ints.
    """
    if not (
        (total_shards is None and shard is None)
        or (isinstance(total_shards, int) and isinstance(shard, int))
    ):
        raise ValueError(f"{shard=} and {total_shards=} must be both None or both ints")

    task_type = Task.IS2RE

    out_dir.mkdir(exist_ok=True, parents=True)

    model_name_ = model_name.split("/")[-1]
    model_name_ = model_name_.replace(":", "-")
    out_path = f"{out_dir}/{model_name_}-{today}.json.gz"

    if total_shards is not None:
        out_path = f"{out_dir}/{model_name_}-{today}-shard-{shard:>03}.json.gz"

    if os.path.isfile(out_path):
        raise SystemExit(f"{out_path=} already exists, exciting early")

    # This is inside the script because accessing the variables causes a download
    # to be triggered if they are not present, meaning it's better to only load them
    # if the script is actually going to be run.
    from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm

    DATA_PATHS = {
        Task.RS2RE: DataFiles.wbm_relaxed_atoms.path,
        Task.IS2RE: DataFiles.wbm_initial_atoms.path,
    }

    data_path = DATA_PATHS[task_type]

    print(f"Loading model {model_name} on {device}")
    model = ORB_PRETRAINED_MODELS[model_name]()
    model.to(device)
    orb_calc = ORBCalculator(model, device=device)

    df_in = pd.read_json(data_path).set_index(str(Key.mat_id))
    if total_shards is not None and shard is not None:
        df_in = np.array_split(df_in, total_shards)[shard - 1]

    run_params = {
        "data_path": data_path,
        Key.task_type: task_type,
        "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
        "max_steps": max_steps,
        "force_max": force_max,
        "ase_optimizer": ase_optimizer,
        "device": device,
        "model_name": model_name,
        "ase_filter": ase_filter,
        "shard": shard,
        "total_shards": total_shards,
    }

    wandb.init(project="matbench-discovery", config=run_params)

    relax_results: dict[str, dict[str, Any]] = {}
    input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]

    if task_type == Task.RS2RE:
        df_in[input_col] = [cse["structure"] for cse in df_in[Key.cse]]

    if limit is not None:
        df_in = df_in.head(limit)

    structs = df_in[input_col].map(Structure.from_dict).to_dict()
    filter_cls = FILTERS[ase_filter]

    for material_id in tqdm(structs, desc="Relaxing"):
        if material_id in relax_results:
            continue
        try:
            atoms = structs[material_id].to_ase_atoms()
            atoms.calc = orb_calc

            if cell_opt:
                atoms = filter_cls(atoms)
            optim_cls = OPTIMIZERS[ase_optimizer]
            optimizer = optim_cls(atoms, logfile="/dev/null")

            optimizer.run(fmax=force_max, steps=max_steps)
            energy = atoms.get_potential_energy()  # relaxed energy
            optimized_structure = AseAtomsAdaptor.get_structure(
                getattr(atoms, "atoms", atoms)  # atoms might be wrapped in ase filter
            )

            relax_results[material_id] = {
                "structure": optimized_structure,
                "energy": energy,
            }

        except Exception as exc:
            print(f"Failed to relax {material_id}: {exc!r}")
            continue

    df_out = pd.DataFrame(relax_results).T.add_prefix("orb_")
    df_out.index.name = str(Key.mat_id)

    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)

    # fix the index name
    df_wbm.index.name = str(df_wbm.index.name)

    df_wbm[PREDICTED_ENERGY_COL] = df_out[PREDICTED_ENERGY_COL]

    predictions_df = df_wbm[
        [str(MbdKey.dft_energy), PREDICTED_ENERGY_COL, str(Key.formula)]
    ]
    predictions_df = predictions_df.reset_index().dropna()

    table = wandb.Table(dataframe=predictions_df)

    title = f"ORB {task_type} ({len(df_out):,})"
    wandb_scatter(
        table,
        fields=dict(x=str(MbdKey.dft_energy), y=PREDICTED_ENERGY_COL),
        title=title,
    )

    wandb.log_artifact(out_path, type=f"orb-wbm-{task_type}")


if __name__ == "__main__":
    app()
