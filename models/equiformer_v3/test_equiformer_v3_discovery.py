"""Perform multiple relaxations in parallel on multiple GPUs.

Example::

    for i in $(seq 0 31); do
        device=$((i % 8))
        CUDA_VISIBLE_DEVICES=$device python test_discovery.py \
            --checkpoint-path $CHECKPOINT_PATH \
            --output-path $OUTPUT_DIR \
            --data-path $DATA_PATH \
            --num-jobs 32 --job-index $i &
    done

To terminate all parallel relaxations::

    pkill -f "python.*test_discovery.py"
"""

import json
import os
import random
from importlib.metadata import version
from pathlib import Path
from typing import Annotated, Any, Literal

import numpy as np
import pandas as pd
import torch
import tqdm
import typer
from ase import Atoms
from ase.filters import Filter, FrechetCellFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from fairchem.core import OCPCalculator
from fairchem.core.common.distutils import is_master
from fairchem.core.common.utils import setup_logging
from fairchem.core.datasets import AseDBDataset
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from submitit.helpers import Checkpointable
from torch.utils.data import Subset

from matbench_discovery.data import as_dict_handler, df_wbm

# BASE_PATH = Path("/data/NFS/radish/ylliao/omat24_dataset/matbench_discovery/")
# DATABASE_PATH = {
#    #Task.RS2RE: str(BASE_PATH / "WBM_RS2RE.aselmdb"),
#    Task.IS2RE: str(BASE_PATH / "WBM_IS2RE.aselmdb"),
# }

FILTER_CLS: dict[str, type[Filter]] = {
    "frechet": FrechetCellFilter,
    "unit": UnitCellFilter,
}
OPTIM_CLS: dict[str, type[Optimizer]] = {"FIRE": FIRE, "LBFGS": LBFGS, "BFGS": BFGS}


_rng = np.random.default_rng(0)


def seed_everywhere(seed: int) -> None:
    global _rng  # noqa: PLW0603
    random.seed(seed)
    _rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def setup_distributed() -> tuple[int, int, int] | None:
    if "RANK" in os.environ:
        rank = int(os.environ["RANK"])
        world_size = int(os.environ["WORLD_SIZE"])
        local_rank = int(os.environ.get("LOCAL_RANK", rank))
        if not torch.distributed.is_initialized():
            torch.distributed.init_process_group(
                backend="nccl",
                rank=rank,
                world_size=world_size,
            )
        return rank, world_size, local_rank
    return None


class AseDBSubset(Subset):
    def get_atoms(self, idx: int) -> Atoms:
        return self.dataset.get_atoms(self.indices[idx])


class RelaxJob(Checkpointable):
    """Submitit checkpointable MLFF relax job to handle preemptions gracefully"""

    def __init__(self) -> None:
        self.relax_results: dict[str, dict[str, Any]] = {}

    def __call__(
        self,
        data_path: Path,
        checkpoint_path: Path,
        output_path: Path,
        optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE",
        cell_filter: Literal["frechet", "unit"] | None = "frechet",
        force_max: float = 0.02,
        max_steps: int = 500,
        optimizer_params: dict[str, Any] | None = None,
        device: Literal["cuda", "cpu"] = "cuda",
        *,
        use_amp: bool = False,
        num_jobs: int = 1,
        job_index: int = 0,
        debug: bool = False,
    ) -> None:
        self.output_path = output_path
        self.num_jobs = num_jobs
        self.job_index = job_index

        seed_everywhere(0)

        run_name = f"index@{self.job_index}-total@{self.num_jobs}"
        json_path = output_path / ("results_" + run_name + ".json.gz")

        print(f"Starting ASE relaxation job {run_name}.")

        calc = OCPCalculator(
            checkpoint_path=checkpoint_path, cpu=(device == "cpu"), seed=0
        )
        if use_amp is False:  # disable the scaler
            calc.trainer.scaler = None

        dataset = AseDBDataset(dict(src=str(data_path)))

        if debug:
            indices = _rng.permutation(len(dataset))
            indices = indices[0:1000]
            indices = np.array_split(indices, num_jobs)[job_index]
        else:
            # indices = np.array_split(range(len(dataset)), num_jobs)[job_index]
            temp = list(range(len(dataset)))
            indices = [temp[i::num_jobs] for i in range(num_jobs)]
            indices = indices[job_index]
        dataset = AseDBSubset(dataset, indices)

        optimizer_params = optimizer_params or {}
        run_params = {
            "data_path": str(data_path),
            "versions": {
                dep: version(dep) for dep in ("fairchem-core", "numpy", "torch")
            },
            Key.task_type: "IS2RE",
            "max_steps": max_steps,
            "force_max": force_max,
            "device": device,
            Key.model_params: sum(p.numel() for p in calc.trainer.model.parameters()),
            "optimizer": optimizer,
            "filter": cell_filter,
            "optimizer_params": optimizer_params,
        }
        with open(output_path / "run_params.json", mode="w") as file:
            json.dump(run_params, file, indent=4)

        self._ase_relax(
            dataset=dataset,
            calculator=calc,
            optimizer=optimizer,
            cell_filter=cell_filter,
            force_max=force_max,
            max_steps=max_steps,
            optimizer_params=optimizer_params,
        )

        df_out = pd.DataFrame(self.relax_results).T.add_prefix("")
        df_out.index.name = Key.mat_id
        df_out.reset_index().to_json(
            json_path, default_handler=as_dict_handler, orient="records", lines=True
        )
        e_pred_col = "energy"
        df_wbm[e_pred_col] = df_out[e_pred_col]
        return

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

        for idx in tqdm.tqdm(
            range(len(dataset)),
            position=self.job_index,
            desc=f"job_index {self.job_index}",
        ):
            atoms = dataset.get_atoms(idx)
            material_id = atoms.info["material_id"]
            if material_id in self.relax_results:
                print(f"Structure {material_id} has already been relaxed.")
                continue
            try:
                atoms.calc = calculator

                filtered_atoms = atoms if filter_cls is None else filter_cls(atoms)
                optim_inst = optim_cls(
                    filtered_atoms,
                    logfile="/dev/null",
                    **optimizer_params,
                )

                optim_inst.run(fmax=force_max, steps=max_steps)

                energy = atoms.get_potential_energy()
                unwrapped = getattr(filtered_atoms, "atoms", atoms)
                structure = AseAtomsAdaptor.get_structure(unwrapped)

                self.relax_results[material_id] = {
                    "structure": structure,
                    "energy": energy,
                }
            except Exception:
                print(f"Failed to relax {material_id}")
                output_err_filename = (
                    f"error_index@{self.job_index}_total@{self.num_jobs}.txt"
                )
                with open((self.output_path / output_err_filename), "a") as txt_file:
                    txt_file.write(f"{material_id}\n")
                continue


def run_relax(
    data_path: Annotated[Path, typer.Option()],
    checkpoint_path: Annotated[Path, typer.Option()],
    output_path: Annotated[
        Path, typer.Option(help="Output path to write results files")
    ],
    optimizer: Annotated[
        Literal["FIRE", "LBFGS", "BFGS"],
        typer.Option(help="Optimizer for relaxations: 'FIRE', 'BFGS', 'LBFGS'"),
    ] = "FIRE",
    cell_filter: Annotated[
        Literal["frechet", "unit"] | None,
        typer.Option(help="Filter for cell relaxation"),
    ] = "frechet",
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
    device: Annotated[
        Literal["cuda", "cpu"], typer.Option(help="Device to use torch")
    ] = "cuda",
    *,
    use_amp: Annotated[
        bool, typer.Option(help="Use automatic mixed precision")
    ] = False,
    num_jobs: Annotated[
        int, typer.Option(help="Number of parallel jobs to relax structures")
    ] = 1,
    job_index: Annotated[int, typer.Option(help="Index of jobs")] = 0,
    debug: Annotated[
        bool, typer.Option(help="Debug mode and evaluate on a subset of structures")
    ] = False,
) -> None:
    setup_logging()
    # rank, world_size, local_rank = setup_distributed()
    # torch.cuda.set_device(local_rank)

    if is_master():
        os.makedirs(output_path, exist_ok=True)

    RelaxJob()(
        data_path,
        checkpoint_path,
        output_path,
        optimizer,
        cell_filter,
        force_max,
        max_steps,
        {} if optimizer_params is None else json.loads(optimizer_params),
        device,
        use_amp=use_amp,
        num_jobs=num_jobs,
        job_index=job_index,
        debug=debug,
    )
    return


if __name__ == "__main__":
    typer.run(run_relax)
