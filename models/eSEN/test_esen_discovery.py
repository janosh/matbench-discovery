"""
Copyright (c) Meta Platforms, Inc. and its affiliates.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

from __future__ import annotations

import random
from importlib.metadata import version
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import torch
import wandb
from ase.filters import FrechetCellFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, LBFGS
from fairchem.core import OCPCalculator
from fairchem.core.datasets import AseDBDataset
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm, trange

from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.enums import MbdKey


def seed_everywhere(seed: int) -> None:
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


BASE_PATH = Path("Matbench-Discovery_data_rootpath")

DATABASE_PATH = {
    "is2re": str(BASE_PATH / "WBM_IS2RE.aselmdb"),
}

FILTER_CLS = {"frechet": FrechetCellFilter, "unit": UnitCellFilter}
OPTIM_CLS = {"FIRE": FIRE, "LBFGS": LBFGS, "BFGS": BFGS}


class MBDRunner:
    """
    we will create a separate script for running MACE.
    """

    def __init__(
        self,
        seed: int,
        model_dir: str,
        save_name: str,
        identifier: str,
        task_type: str = "is2re",
        optimizer: Literal["FIRE", "LBFGS", "BFGS"] = "FIRE",
        cell_filter: Literal["frechet", "unit"] = "frechet",
        force_max: float = 0.02,
        max_steps: int = 500,
        optimizer_params: dict[str, Any] | None = None,
        device: Literal["cuda", "cpu"] = "cuda",
        batch_size: int = 1,
        num_jobs: int = 1,
        num_workers: int = 0,
    ) -> None:
        self.seed = seed
        self.model_dir = model_dir
        self.save_name = save_name
        self.identifier = identifier

        self.task_type = task_type
        self.optimizer = optimizer
        self.cell_filter = cell_filter
        self.force_max = force_max
        self.max_steps = max_steps
        self.optimizer_params = optimizer_params
        self.device = device
        self.batch_size = batch_size

        self.num_jobs = num_jobs
        self.num_workers = num_workers

    def run(self, job_number: int = 0) -> None:
        self.relax_results: dict[str, Any] = {}

        save_dir = (
            Path(self.model_dir) / self.save_name / f"{self.identifier}_{self.seed}"
        )
        (save_dir).mkdir(parents=True, exist_ok=True)
        self.save_dir = save_dir

        seed_everywhere(self.seed)
        model_ckpt = str(Path(self.model_dir) / "checkpoint.pt")
        calc = OCPCalculator(checkpoint_path=model_ckpt, cpu=False, seed=0)
        calc.trainer.scaler = None
        num_model_params = sum(p.numel() for p in calc.trainer.model.parameters())

        data_path = DATABASE_PATH[self.task_type]
        dataset = AseDBDataset(dict(src=data_path))

        indices = None
        if self.num_jobs > 1:
            indices = np.array_split(range(len(dataset)), self.num_jobs)[job_number]

        optimizer_params = self.optimizer_params or {}
        run_params = {
            "data_path": data_path,
            "versions": {
                dep: version(dep) for dep in ("fairchem-core", "numpy", "torch")
            },
            Key.task_type: self.task_type,
            "max_steps": self.max_steps,
            "force_max": self.force_max,
            "device": self.device,
            Key.model_params: num_model_params,
            "optimizer": self.optimizer,
            "filter": self.cell_filter,
            "optimizer_params": self.optimizer_params,
        }
        wandb.init(
            project="matbench-discovery", config=run_params, name=self.identifier
        )

        self._ase_relax(
            dataset=dataset,
            indices=indices,
            calculator=calc,
            optimizer_cls=self.optimizer,
            cell_filter=self.cell_filter,
            force_max=self.force_max,
            max_steps=self.max_steps,
            optimizer_params=optimizer_params,
        )

        df_out = pd.DataFrame(self.relax_results).T
        df_out.index.name = Key.mat_id
        df_out.reset_index().to_json(
            save_dir / f"{self.identifier}_{job_number}.json.gz",
            default_handler=as_dict_handler,
            orient="records",
            lines=True,
        )
        e_pred_col = "pred_energy"
        df_wbm[e_pred_col] = df_out[e_pred_col]

        # join prediction
        file_paths = list(self.save_dir.glob(f"{self.identifier}_*.json.gz"))
        num_job_finished = len(file_paths)
        if num_job_finished == self.num_jobs:
            self.join_prediction(file_paths)

    def _ase_relax(
        self,
        dataset: AseDBDataset,
        indices: np.ndarray | None,
        calculator: OCPCalculator,
        optimizer_cls: Literal["FIRE", "LBFGS", "BFGS"],
        cell_filter: Literal["frechet", "unit"],
        force_max: float,
        max_steps: int,
        optimizer_params: dict[str, Any],
    ) -> None:
        """Run WBM relaxations using an ASE optimizer."""
        filter_cls = FILTER_CLS.get(cell_filter)
        optim_cls = OPTIM_CLS[optimizer_cls]

        # iterate over indices if provided, otherwise over full dataset
        iteration_indices = indices if indices is not None else range(len(dataset))

        for idx in trange(len(iteration_indices), desc="Relaxing with ASE"):
            dataset_idx = int(iteration_indices[idx])
            atoms = dataset.get_atoms(dataset_idx)
            material_id = atoms.info["sid"]
            if material_id in self.relax_results:
                continue
            try:
                atoms.calc = calculator

                if filter_cls is not None:
                    optimizer = optim_cls(
                        filter_cls(atoms), logfile="/dev/null", **optimizer_params
                    )
                else:
                    optimizer = optim_cls(
                        atoms, logfile="/dev/null", **optimizer_params
                    )

                optimizer.run(fmax=force_max, steps=max_steps)

                energy = atoms.get_potential_energy()
                structure = AseAtomsAdaptor.get_structure(atoms)

                self.relax_results[material_id] = {
                    "pred_structure": structure,
                    "pred_energy": energy,
                }
            except Exception:
                continue

    def join_prediction(self, file_paths: list[Path] | None = None) -> None:
        dfs: dict[str, pd.DataFrame] = {}
        for file_path in tqdm(file_paths, desc="Loading prediction files"):
            if file_path in dfs:
                continue
            dfs[file_path] = pd.read_json(file_path).set_index(Key.mat_id)

        df_fairchem = pd.concat(dfs.values()).round(4)

        # %%
        df_wbm_cse = pd.read_json(
            DataFiles.wbm_computed_structure_entries.path, lines=True
        ).set_index(Key.mat_id)
        df_wbm_cse[Key.computed_structure_entry] = [
            ComputedStructureEntry.from_dict(dct)
            for dct in tqdm(
                df_wbm_cse[Key.computed_structure_entry], desc="Creating pmg CSEs"
            )
        ]

        # corrections applied below are structure-dependent (for oxides and sulfides)
        cse: ComputedStructureEntry
        for row in tqdm(
            df_fairchem.itertuples(), total=len(df_fairchem), desc="ML energies to CSEs"
        ):
            mat_id, struct_dict, pred_energy, *_ = row
            mlip_struct = Structure.from_dict(struct_dict)
            cse = df_wbm_cse.loc[mat_id, Key.computed_structure_entry]
            cse._energy = pred_energy  # noqa: SLF001
            cse._structure = mlip_struct  # noqa: SLF001
            df_fairchem.loc[mat_id, Key.computed_structure_entry] = cse

        # %% apply energy corrections
        processed = MaterialsProject2020Compatibility().process_entries(
            df_fairchem[Key.computed_structure_entry], verbose=True, clean=True
        )
        if len(processed) != len(df_fairchem):
            raise ValueError(
                f"not all entries processed: {len(processed)=} {len(df_fairchem)=}"
            )

        # %% compute corrected formation energies
        df_fairchem[Key.formula] = df_wbm[Key.formula]
        df_fairchem["pred_e_form_per_atom"] = [
            get_e_form_per_atom(dict(energy=cse.energy, composition=formula))
            for formula, cse in tqdm(
                df_fairchem.set_index(Key.formula)[
                    Key.computed_structure_entry
                ].items(),
                total=len(df_fairchem),
                desc="Computing formation energies",
            )
        ]
        df_wbm[[*df_fairchem]] = df_fairchem

        # %%
        bad_mask = abs(df_wbm["pred_e_form_per_atom"] - df_wbm[MbdKey.e_form_dft]) > 5
        n_preds = len(df_wbm["pred_e_form_per_atom"].dropna())
        print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
        df_fairchem = df_fairchem.round(4)

        df_fairchem.select_dtypes("number").to_csv(
            self.save_dir / f"{self.identifier}.csv.gz"
        )

        df_bad = df_fairchem[bad_mask].drop(
            columns=[Key.computed_structure_entry, "pred_structure"]
        )
        df_bad[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
        df_bad.to_csv(self.save_dir / "bad.csv")

        df_fairchem.reset_index().to_json(
            self.save_dir / f"{self.identifier}.json.gz",
            default_handler=as_dict_handler,
            orient="records",
            lines=True,
        )
