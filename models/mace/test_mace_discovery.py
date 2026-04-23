# %%
import os
from importlib.metadata import version
from typing import Any, Final

import pandas as pd
import plotly.express as px
import torch
import wandb
from ferrox import SimLogLevel, Structure
from ferrox.mlff import MaceModel
from ferrox.relax import relax
from ferrox.trajectory import Trajectory
from mace.calculators.foundations_models import download_mace_mp_checkpoint
from mace.tools import count_parameters
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import hpc, timestamp, today
from matbench_discovery.data import as_dict_handler, df_wbm, structures_from_zip
from matbench_discovery.enums import DataFiles, MbdKey, Model, Task
from matbench_discovery.plots import wandb_scatter

__author__ = "Janosh Riebesell"
__date__ = "2024-12-09"


# %%
smoke_test = False
if smoke_test:
    print(f"Warning: {smoke_test=}, will not write relaxed structures to disk!")
task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 52
optimizer_name = "fire"
device = "cuda" if torch.cuda.is_available() else "cpu"
# whether to record intermediate structures into ferrox Trajectory
record_traj = False  # skip trajectory recording for full benchmark run
model_name = os.getenv("MODEL_NAME", Model.mace_mpa_0)
job_name = f"{model_name}/{today}-wbm-{task_type}-{optimizer_name.upper()}"
out_dir = f"{module_dir}/{job_name}"
os.makedirs(out_dir, exist_ok=True)
checkpoint_urls: Final[set[str]] = {
    "https://github.com/ACEsuit/mace-foundations/releases/download/mace_omat_0/mace-omat-0-medium.model",
    "https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model",
    "https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/mace-mpa-0-medium.model",
    "https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0/2023-12-03-mace-128-L1_epoch-199.model",
}
checkpoint = {url.split("/")[-1].rsplit(".model")[0]: url for url in checkpoint_urls}[
    model_name
]
print(f"{model_name=}")

slurm_vars = hpc.slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    array=f"1-{slurm_array_task_count}",
    slurm_flags="--gpus=1 --partition=h100 --cpus-per-task=4 --mem=32G --time=02:00:00",
    pre_cmd="source ~/periodic-mono/.venv/bin/activate && source ~/.cargo/env",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "1"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exiting early")


# %%
data_path = {
    Task.RS2RE: DataFiles.wbm_relaxed_atoms.path,
    Task.IS2RE: DataFiles.wbm_initial_atoms.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
e_pred_col = "mace_energy"
max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A
dtype = "float64"

# Download checkpoint and load model via ferrox (no ASE calculator needed)
model_path = download_mace_mp_checkpoint(checkpoint)
raw_model = torch.load(model_path, map_location="cpu", weights_only=False)
if dtype == "float64":
    raw_model = raw_model.double()
raw_model = raw_model.to(device)
raw_model.eval()
mace_model = MaceModel.from_model(raw_model, device)

print(f"Read data from {data_path}")
struct_dict = structures_from_zip(data_path)

if slurm_array_job_id == "debug":
    if smoke_test:
        struct_dict = dict(list(struct_dict.items())[:128])
elif slurm_array_task_count > 1:
    structures = list(struct_dict.values())
    chunks = hpc.chunk_by_lens(structures, n_chunks=slurm_array_task_count)
    my_chunk_ids = {id(struct) for struct in chunks[slurm_array_task_id - 1]}
    struct_dict = {
        mat_id: struct
        for mat_id, struct in struct_dict.items()
        if id(struct) in my_chunk_ids
    }


# %%
run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("mace-torch", "numpy", "torch")},
    "checkpoint": checkpoint,
    Key.task_type: task_type,
    "n_structures": len(struct_dict),
    "slurm_vars": slurm_vars,
    "max_steps": max_steps,
    "record_traj": record_traj,
    "force_max": force_max,
    "optimizer": optimizer_name,
    "device": device,
    Key.model_params: count_parameters(mace_model.model),
    "model_name": model_name,
    "dtype": dtype,
    "cell_filter": "FrechetCellFilter",
}

run_name = f"{job_name}-{slurm_array_task_id}"

wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %% time
relax_results: dict[str, dict[str, Any]] = {}

for mat_id, structure in tqdm(struct_dict.items(), desc="Relaxing"):
    if mat_id in relax_results:
        continue
    try:
        result = relax(
            structure,
            mace_model,
            optimizer=optimizer_name,
            fmax=force_max,
            max_steps=max_steps,
            relax_cell=True,
            log=SimLogLevel.NONE,
        )
        relaxed_struct: Structure = result.final_structure
        energy: float = result.final_energy
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}

        if record_traj and result.history:
            # Build ferrox Trajectory from optimization history
            coords = [step.positions.tolist() for step in result.history]
            lattice_rows: list[list[float]] = []
            for step in result.history:
                if step.cell is not None:
                    lattice_rows.extend(step.cell.tolist())
            energies_list = [step.energy for step in result.history]

            mace_traj = Trajectory(
                species=structure.species_strings,
                coords=coords,
                lattice=lattice_rows,
                constant_lattice=False,
                frame_properties=[{"energy": eng} for eng in energies_list],
            )
            relax_results[mat_id]["trajectory"] = mace_traj
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue


# %%
df_out = pd.DataFrame(relax_results).T.add_prefix("mace_")
df_out.index.name = Key.mat_id
if not smoke_test:
    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )


# %%
if Key.trajectory in df_out:
    energy_series = df_out[Key.trajectory].map(
        lambda x: [d["energy"] / len(x.species) for d in x.frame_properties]
    )

    # Create a DataFrame from the Series
    df_energies = pd.DataFrame(energy_series.tolist()).T
    df_energies.columns = df_out.index
    df_energies["Step"] = df_energies.index

    # Melt the DataFrame to long format
    df_melted = df_energies.melt(
        id_vars=["Step"], var_name="Trajectory", value_name="Energy"
    )

    # Create the line plot
    fig = px.line(
        df_melted,
        x="Step",
        y="Energy",
        color="Trajectory",
        title="Trajectory Energies",
        labels={"Step": "Optimization Step", "Energy": "Energy"},
        line_group="Trajectory",
    )

    # Customize the layout if needed
    fig.update_layout(
        xaxis_title="Optimization Step", yaxis_title="Energy", legend_title="Trajectory"
    )

    # Show the plot
    fig.show()


# %%
df_wbm[e_pred_col] = df_out[e_pred_col]

table = wandb.Table(
    dataframe=df_wbm[[MbdKey.dft_energy, e_pred_col, Key.formula]]
    .reset_index()
    .dropna()
)

title = f"MACE {task_type} ({len(df_out):,})"
wandb_scatter(table, fields=dict(x=MbdKey.dft_energy, y=e_pred_col), title=title)

wandb.log_artifact(out_path, type=f"mace-wbm-{task_type}")
