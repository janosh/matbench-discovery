import random
from importlib.metadata import version
from pathlib import Path
from typing import Any, Literal
import numpy as np
import pandas as pd
import swanlab
import torch
from ase import Atoms
from ase.filters import FrechetCellFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, LBFGS
from matnova.core.common.relaxation.ase_utils import OCPCalculator
from matnova.core.datasets import AseDBDataset
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from torch.utils.data import Subset
from tqdm import tqdm, trange
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.enums import MbdKey

def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None

class AseDBSubset(Subset):
    def get_atoms(self, idx: int) -> Atoms:
        return self.dataset.get_atoms(self.indices[idx])
    
def seed_everywhere(seed: int) -> None:
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

ROOT = "/mnt/public/LiuLiu/fairchem/dataset/wbm"
DATABASE_PATH = {"is2re": f"{ROOT}/WBM_IS2RE.lmdb"}
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

        save_dir = Path(self.model_dir) / f"{self.identifier}_{self.seed}"
        (save_dir).mkdir(parents=True, exist_ok=True)
        self.save_dir = save_dir

        seed_everywhere(self.seed)
        model_ckpt = str(Path(self.model_dir) / "checkpoint.pt")
        calc = OCPCalculator(checkpoint_path=model_ckpt, cpu=False, seed=0)

        # calc.trainer.scaler = None
        num_model_params = sum(p.numel() for p in calc.trainer.model.parameters())

        data_path = DATABASE_PATH[self.task_type]
        dataset = AseDBDataset(dict(src=data_path))


        if self.num_jobs > 1:
            indices = np.array_split(range(min(100, len(dataset))), self.num_jobs)[
                job_number
           ]
            dataset = AseDBSubset(dataset, indices)

        optimizer_params = self.optimizer_params or {}
        run_params = {
            "data_path": data_path,
            "versions": {dep: version(dep) for dep in ("numpy", "torch")},
            Key.task_type: self.task_type,
            "max_steps": self.max_steps,
            "force_max": self.force_max,
            "device": self.device,
            Key.model_params: num_model_params,
            "optimizer": self.optimizer,
            "filter": self.cell_filter,
            "optimizer_params": self.optimizer_params,
        }

        swanlab.init(project="matbench-discovery", config=run_params)

        self._ase_relax(
            dataset=dataset,
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
        
        self._save_relaxed_structures_to_csv(df_out)

        df_wbm[e_pred_col] = df_out[e_pred_col]

        file_paths = list(self.save_dir.glob(f"{self.identifier}_*.json.gz"))
        num_job_finished = len(file_paths)

        if num_job_finished == self.num_jobs:
            self.join_prediction(file_paths)


    def _ase_relax(
        self,
        dataset: AseDBDataset | AseDBSubset,
        calculator: OCPCalculator,
        optimizer_cls: Literal["FIRE", "LBFGS", "BFGS"],
        cell_filter: Literal["frechet", "unit"],
        force_max: float,
        max_steps: int,
        optimizer_params: dict[str, Any],
    ) -> None:
        """Run WBM relaxations using an ASE optimizer."""
        filter_cls = FILTER_CLS.get(cell_filter)
        optim_cls = OPTIM_CLS[optimizer_cls]  # 'ase.optimize.fire.FIRE'

        skip_num = 0
        for i in trange(len(dataset), desc="Relaxing with ASE"):
            
            atoms = dataset.get_atoms(i)
            material_id = atoms.info["material_id"]
            if material_id in self.relax_results:
                skip_num += 1
                print(f"Skipping {material_id}, already relaxed")
                continue
            try:
                atoms.calc = calculator

                if filter_cls is not None:
                    optimizer = optim_cls(
                        filter_cls(atoms), 
                        logfile="/dev/null", 
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
            except Exception as e:
                skip_num += 1
                print(f"Failed to relax {material_id}: {e}")
                continue

    def join_prediction(self, file_paths: list[Path] | None = None) -> None:
        dfs: dict[str, pd.DataFrame] = {}
        
        for file_path in tqdm(file_paths, desc="Loading prediction files"):
            if file_path in dfs:
                continue
            read_file = pd.read_json(file_path, lines=True)
            dfs[file_path] = read_file.set_index(Key.mat_id.value)

        df_fairchem = pd.concat(dfs.values()).round(4)
        self._save_relaxed_structures_to_csv(df_fairchem)
        
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

        cse: ComputedStructureEntry
        for row in tqdm(
            df_fairchem.itertuples(), total=len(df_fairchem), desc="ML energies to CSEs"
        ):
            *_, mat_id, struct_dict, pred_energy = row
            mlip_struct = Structure.from_dict(struct_dict)
            cse = df_wbm_cse.loc[mat_id, Key.computed_structure_entry]
            cse._energy = pred_energy  # noqa: SLF001
            cse._structure = mlip_struct  # noqa: SLF001
            df_fairchem.loc[mat_id, Key.computed_structure_entry] = cse

        entries_col = df_fairchem[Key.computed_structure_entry]

         mask = entries_col.notna() & entries_col.apply(
            lambda x: isinstance(x, ComputedStructureEntry)
         )

        df_fairchem = df_fairchem[mask]

        processed = MaterialsProject2020Compatibility().process_entries(
            df_fairchem[Key.computed_structure_entry], verbose=True, clean=True
        )

        if len(processed) != len(df_fairchem):
            raise ValueError(
                f"not all entries processed: {len(processed)=} {len(df_fairchem)=}"
            )

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
        temp_mask = abs(df_wbm["pred_e_form_per_atom"] - df_wbm[MbdKey.e_form_dft])
        bad_mask = temp_mask > 5
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

    def _save_relaxed_structures_to_csv(self, df_out: pd.DataFrame) -> None:
        try:
            csv_data = []
            for material_id, row in df_out.iterrows():
                if "pred_structure" in row and "pred_energy" in row:
                    energy = row["pred_energy"]
                    structure_info = {
                        "material_id": material_id,
                        "pred_energy": energy,
                    }

                    csv_data.append(structure_info)

            if csv_data:
                df_csv = pd.DataFrame(csv_data)
                csv_filename = f"{self.identifier}_relaxed_structures.csv"
                df_csv.to_csv(csv_filename, index=False)

        except Exception as e:
            print(f"error: {e}")
