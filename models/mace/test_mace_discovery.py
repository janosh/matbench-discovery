# %%
import os
from copy import deepcopy
from importlib.metadata import version
from typing import Any, Final

import pandas as pd
import plotly.express as px
import wandb
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from mace.calculators import mace_mp
from mace.tools import count_parameters
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import hpc, timestamp, today
from matbench_discovery.data import as_dict_handler, ase_atoms_from_zip, df_wbm
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
slurm_array_task_count = 200
ase_optimizer = "FIRE"
# device = "cuda" if torch.cuda.is_available() else "cpu"
device = "cpu"
# whether to record intermediate structures into pymatgen Trajectory
record_traj = True  # has no effect if relax_cell is False
model_name = os.getenv("MODEL_NAME", Model.mace_mp_0)
job_name = f"{model_name}/{today}-wbm-{task_type}-{ase_optimizer}"
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
    # slurm_flags="--qos shared --constraint gpu --gpus 1",
    slurm_flags="--ntasks=1 --cpus-per-task=1 --partition high-priority",
    submit_as_temp_file=False,
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
mace_calc = mace_mp(model=checkpoint, device=device, default_dtype=dtype)

print(f"Read data from {data_path}")
atoms_list = ase_atoms_from_zip(data_path)

if slurm_array_job_id == "debug":
    if smoke_test:
        atoms_list = atoms_list[:128]
    else:
        pass
elif slurm_array_task_count > 1:
    atoms_list = hpc.chunk_by_lens(atoms_list, n_chunks=slurm_array_task_count)[
        slurm_array_task_id - 1
    ]


# %%
run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("mace-torch", "numpy", "torch")},
    "checkpoint": checkpoint,
    Key.task_type: task_type,
    "n_structures": len(atoms_list),
    "slurm_vars": slurm_vars,
    "max_steps": max_steps,
    "record_traj": record_traj,
    "force_max": force_max,
    "ase_optimizer": ase_optimizer,
    "device": device,
    Key.model_params: count_parameters(mace_calc.models[0]),
    "model_name": model_name,
    "dtype": dtype,
    "cell_filter": "FrechetCellFilter",
}

run_name = f"{job_name}-{slurm_array_task_id}"

wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %% time
relax_results: dict[str, dict[str, Any]] = {}
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

for atoms in tqdm(deepcopy(atoms_list), desc="Relaxing"):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = mace_calc
        if max_steps > 0:
            filtered_atoms = FrechetCellFilter(atoms)
            optimizer = optim_cls(filtered_atoms, logfile=None)

            if record_traj:
                coords, lattices, energies = [], [], []
                # attach observer functions to the optimizer
                optimizer.attach(lambda: coords.append(atoms.get_positions()))  # noqa: B023
                optimizer.attach(lambda: lattices.append(atoms.get_cell()))  # noqa: B023
                optimizer.attach(lambda: energies.append(atoms.get_potential_energy()))  # noqa: B023

            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by FrechetCellFilter, so need to getattr
        relaxed_struct = AseAtomsAdaptor.get_structure(atoms)
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}

        coords = locals().get("coords", [])
        lattices = locals().get("lattices", [])
        energies = locals().get("energies", [])
        if record_traj and coords and lattices and energies:
            mace_traj = Trajectory(
                species=atoms.get_chemical_symbols(),
                coords=coords,
                lattice=lattices,
                constant_lattice=False,
                frame_properties=[{"energy": energy} for energy in energies],
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
