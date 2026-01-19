# %%
import sys

sys.path.insert(0, "/data/fywang/code/fairchem-tracegrad-hanhai/src")
import json
import logging
import os
from pathlib import Path
from typing import Annotated, Any, Literal

import numpy as np
import pandas as pd
import torch
import typer
import wandb
from ase import Atoms
from ase.filters import Filter, FrechetCellFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from fairchem.core import OCPCalculator
from fairchem.core.common.utils import setup_logging
from fairchem.core.datasets import AseDBDataset
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from torch.utils.data import Subset
from tqdm import trange

from matbench_discovery import timestamp, today
from matbench_discovery.data import as_dict_handler, df_wbm
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter

BASE_PATH = Path("/data/fywang/code/matbench-discovery/data/wbm")
DATABASE_PATH = {
    # Task.RS2RE: str(BASE_PATH / "WBM_RS2RE.aselmdb"),
    Task.IS2RE: str(BASE_PATH / "WBM_IS2RE_initial_atoms.aselmdb"),
}


# geo_fairchem = pd.read_json("/data/fywang/code/matbench-discovery/models/eqV2/eqV2-s-dens-mp_2024-10-18-wbm-geo-opt.json.gz", compression='gzip')
# discovery_fairchem = pd.read_csv("/data/fywang/code/matbench-discovery/models/eqV2/eqV2-s-dens-mp.csv.gz", compression='gzip' )

FILTER_CLS: dict[str, Filter] = {"frechet": FrechetCellFilter, "unit": UnitCellFilter}
OPTIM_CLS: dict[str, Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS, "BFGS": BFGS}


class AseDBSubset(Subset):
    def get_atoms(self, idx: int) -> Atoms:
        return self.dataset.get_atoms(self.indices[idx])


class RelaxJob:
    """Local execution of the MLFF relax job"""

    def __init__(self) -> None:
        self.relax_results: dict[str, dict[str, Any]] = {}

    def __call__(
        self,
        job_number: int,
        checkpoint_path: Path,
        out_path: Path,
        model_name: str,
        job_name: str,
        *,
        task_type: Task = Task.IS2RE,
        optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE",
        cell_filter: Literal["frechet", "unit"] | None = "frechet",
        force_max: float = 0.02,
        max_steps: int = 500,
        optimizer_params: dict[str, Any] | None = None,
        device: Literal["cuda", "cpu"] = "cuda",
        use_amp: bool = True,
        num_jobs: int = 1,
        debug: bool = False,
    ) -> None:
        run_name = f"{job_name}-{job_number}"
        out_path = out_path / f"{model_name}-{today}-{job_number:>03}.json.gz"

        logging.info(f"Starting ASE relaxation job {run_name}.")

        calc = OCPCalculator(
            checkpoint_path=checkpoint_path,
            cpu=device == "cpu",
            seed=0,
            # trainer="equiformerv2_dens",
        )

        if not use_amp:  # disable the scaler
            calc.trainer.scaler = None

        data_path = DATABASE_PATH[task_type]
        dataset = AseDBDataset(dict(src=data_path))

        if debug:
            indices = np.array_split(range(len(dataset)), 100000)[job_number]
            dataset = AseDBSubset(dataset, indices)
        elif num_jobs > 1:
            indices = np.array_split(range(len(dataset)), num_jobs)[job_number]
            dataset = AseDBSubset(dataset, indices)

        optimizer_params = optimizer_params or {}
        run_params = {
            "data_path": data_path,
            # "versions": {
            #     dep: version(dep) for dep in ("fairchem-core", "numpy", "torch")
            # },
            "task_type": task_type,
            "max_steps": max_steps,
            "force_max": force_max,
            "device": device,
            "model_params": sum(p.numel() for p in calc.trainer.model.parameters()),
            "model_name": model_name,
            "optimizer": optimizer,
            "filter": cell_filter,
            "optimizer_params": optimizer_params,
        }
        wandb.init(project="matbench-discovery", config=run_params, name=run_name)

        self._ase_relax(
            dataset=dataset,
            calculator=calc,
            optimizer=optimizer,
            cell_filter=cell_filter,
            force_max=force_max,
            max_steps=max_steps,
            optimizer_params=optimizer_params,
        )

        df_out = pd.DataFrame(self.relax_results).T.add_prefix(f"{model_name}_")
        df_out.index.name = "material_id"
        df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)
        e_pred_col = f"{model_name}_energy"
        df_wbm[e_pred_col] = df_out[e_pred_col]
        table = wandb.Table(
            # dataframe=df_wbm[["dft_energy", e_pred_col, "formula"]]
            dataframe=df_wbm[[MbdKey.dft_energy, e_pred_col, Key.formula]]
            .reset_index()
            .dropna()
        )

        title = f"FAIRChem {model_name} {task_type} ({len(df_out):,})"
        wandb_scatter(table, fields=dict(x="dft_energy", y=e_pred_col), title=title)
        wandb.log_artifact(out_path, type=f"fairchem-{model_name}-wbm-{task_type}")

    def _ase_relax(
        self,
        dataset: AseDBDataset | AseDBSubset,
        calculator: OCPCalculator,
        optimizer: Literal["FIRE", "LBFGS", "BFGS"],
        cell_filter: Literal["frechet", "unit"] | None,
        force_max: float,
        max_steps: int,
        optimizer_params: dict[str, Any],
    ) -> None:
        """Run WBM relaxations using an ASE optimizer."""
        filter_cls = FILTER_CLS.get(cell_filter or "")
        optim_cls = OPTIM_CLS[optimizer]

        for idx in trange(len(dataset), desc="Relaxing with ASE"):
            atoms = dataset.get_atoms(idx)
            # material_id = atoms.info["sid"]
            material_id = atoms.info["material_id"]
            if material_id in self.relax_results:
                logging.info(f"Structure {material_id} has already been relaxed.")
                continue
            try:
                atoms.calc = calculator

                if filter_cls is not None:
                    atoms = filter_cls(atoms)

                optim_inst = optim_cls(atoms, logfile="/dev/null", **optimizer_params)

                optim_inst.run(fmax=force_max, steps=max_steps)

                # geo_fairchem_idx = geo_fairchem.loc[geo_fairchem['material_id'] == material_id]
                # lattice_fairchem = geo_fairchem_idx["eqV2-31M-dens-MP-p5_structure"][idx]["lattice"]["matrix"]
                # lattice_fairchem = np.array(lattice_fairchem)
                # atoms.set_cell(lattice_fairchem)
                if filter_cls is not None:
                    atoms = atoms.atoms
                energy = atoms.get_potential_energy()
                structure = AseAtomsAdaptor.get_structure(atoms)

                self.relax_results[material_id] = {
                    "structure": structure,
                    "energy": energy,
                }
            except Exception:
                logging.exception(f"Failed to relax {material_id}")
                continue


def run_relax(
    checkpoint_path: Annotated[Path, typer.Option()],
    out_path: Annotated[Path, typer.Option(help="Output path to write results files")],
    model_name: Annotated[str, typer.Option(help="Name given to model")],
    *,
    optimizer: Annotated[
        str, typer.Option(help="Optimizer for relaxations: 'FIRE', 'BFGS', 'LBFGS'")
    ] = "LBFGS",
    cell_filter: Annotated[
        str | None, typer.Option(help="Filter for cell relaxation")
    ] = None,
    force_max: Annotated[
        float, typer.Option(help="Force relaxation convergence threshold")
    ] = 0.02,
    max_steps: Annotated[
        int, typer.Option(help="max number of relaxation steps")
    ] = 500,
    optimizer_params: Annotated[
        str | None,
        typer.Option(help="Optimizer parameters as a json string dictionary"),
    ] = None,
    device: Annotated[str, typer.Option(help="Device to use torch")] = "cuda"
    if torch.cuda.is_available()
    else "cpu",
    use_amp: Annotated[bool, typer.Option(help="Use automatic mixed precision")] = True,
    num_jobs: Annotated[int, typer.Option(help="Number of parallel jobs to run")] = 10,
    debug: Annotated[bool, typer.Option(help="debug mode, run only one job")] = False,
) -> None:
    setup_logging()
    task_type = Task.IS2RE
    job_name = f"{model_name}-wbm-{task_type}-{optimizer}"

    os.makedirs(out_path / model_name, exist_ok=True)

    out_path = out_path / model_name / timestamp
    os.makedirs(out_path)

    optimizer_params = {} if optimizer_params is None else json.loads(optimizer_params)
    args = {
        "checkpoint_path": checkpoint_path,
        "out_path": out_path,
        "model_name": model_name,
        "job_name": job_name,
        "task_type": Task.IS2RE,  # 默认值
        "optimizer": optimizer,  # 传入的优化器
        "cell_filter": cell_filter,  # 默认值为 None
        "force_max": force_max,  # 默认值
        "max_steps": max_steps,  # 默认值
        "optimizer_params": optimizer_params,  # 默认为空字典
        "device": device,  # 默认值
        "use_amp": use_amp,  # 默认值
        "num_jobs": num_jobs,  # 默认值
        "debug": debug,  # 默认值
    }

    job = RelaxJob()
    job(0, **args)  # Running the job locally for a single task


if __name__ == "__main__":
    typer.run(run_relax)
