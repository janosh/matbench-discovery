"""Run DPA-4.0-Pro-MPtrj relaxations for Matbench Discovery.

Example::

    for i in $(seq 0 99); do
        CUDA_VISIBLE_DEVICES=$((i % 8)) python test_dpa4_discovery.py \
            --checkpoint-path /path/to/dpa-4.0-pro-mptrj.pt \
            --data-dir /path/to/wbm-shards \
            --output-dir ./dpa4-relax \
            --num-jobs 100 --job-index $i &
    done
"""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Literal, cast

import pandas as pd
from ase.filters import Filter, FrechetCellFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, LBFGS
from deepmd.calculator import DP
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler
from matbench_discovery.enums import Task

FILTER_CLS: dict[str, type[Filter]] = {
    "frechet": FrechetCellFilter,
    "unit": UnitCellFilter,
}
OPTIM_CLS: dict[str, Any] = {"FIRE": FIRE, "LBFGS": LBFGS, "BFGS": BFGS}


def collect_input_shards(data_dir: Path) -> list[Path]:
    """
    Collect WBM shard JSON files.

    Parameters
    ----------
    data_dir : Path
        Directory containing `wbm_data_*.json` files.

    Returns
    -------
    list[Path]
        Sorted shard paths.

    Raises
    ------
    FileNotFoundError
        If no shard files are found.
    """
    shard_paths = sorted(data_dir.glob("wbm_data_*.json"), key=lambda path: int(path.stem.rsplit("_", 1)[-1]))
    if not shard_paths:
        raise FileNotFoundError(f"No WBM shard files found under {data_dir}")
    return shard_paths


class Relaxer:
    """ASE relaxation wrapper around the DeePMD calculator."""

    def __init__(self, checkpoint_path: Path) -> None:
        """
        Initialize the DeePMD calculator.

        Parameters
        ----------
        checkpoint_path : Path
            Path to the frozen DeePMD checkpoint.
        """
        self.calculator = DP(str(checkpoint_path))
        self.ase_adaptor = AseAtomsAdaptor()

    def relax(
        self,
        atoms: Any,
        optimizer: Literal["FIRE", "LBFGS", "BFGS"],
        cell_filter: Literal["frechet", "unit"] | None,
        force_max: float,
        max_steps: int,
    ) -> tuple[dict[str, Any], float]:
        """
        Relax one structure.

        Parameters
        ----------
        atoms : Any
            ASE atoms or pymatgen structure-like object.
        optimizer : {"FIRE", "LBFGS", "BFGS"}
            ASE optimizer name.
        cell_filter : {"frechet", "unit"} | None
            Cell filter name. If None, only atomic positions are relaxed.
        force_max : float
            Force convergence threshold in eV/A.
        max_steps : int
            Maximum number of optimizer steps.

        Returns
        -------
        tuple[dict[str, Any], float]
            Final structure dictionary and potential energy in eV.
        """
        if isinstance(atoms, (Structure, Molecule)):
            atoms = cast("Any", self.ase_adaptor.get_atoms(atoms))

        atoms.calc = self.calculator
        filter_cls = FILTER_CLS.get(cell_filter or "")
        filtered_atoms = atoms if filter_cls is None else filter_cls(atoms)
        optimizer_inst = OPTIM_CLS[optimizer](filtered_atoms, logfile="/dev/null")
        optimizer_inst.run(fmax=force_max, steps=max_steps)

        energy = atoms.get_potential_energy()
        unwrapped_atoms = getattr(filtered_atoms, "atoms", atoms)
        structure = self.ase_adaptor.get_structure(unwrapped_atoms)
        return structure.as_dict(), energy


def relax_shard(
    input_path: Path,
    output_path: Path,
    relaxer: Relaxer,
    model_name: str,
    optimizer: Literal["FIRE", "LBFGS", "BFGS"],
    cell_filter: Literal["frechet", "unit"] | None,
    force_max: float,
    max_steps: int,
) -> None:
    """
    Relax all structures in one WBM shard.

    Parameters
    ----------
    input_path : Path
        Input WBM shard JSON file.
    output_path : Path
        Output JSON-lines file.
    relaxer : Relaxer
        Initialized DPA4 relaxer.
    model_name : str
        Prefix used for output columns.
    optimizer : {"FIRE", "LBFGS", "BFGS"}
        ASE optimizer name.
    cell_filter : {"frechet", "unit"} | None
        Cell filter name.
    force_max : float
        Force convergence threshold in eV/A.
    max_steps : int
        Maximum number of optimizer steps.
    """
    # === Step 1. Load Structures ===
    input_col = str({Task.IS2RE: Key.initial_struct, Task.RS2RE: Key.final_struct}[Task.IS2RE])
    df_in = pd.read_json(input_path)
    if input_col not in df_in:
        input_col = "initial_structure"
    structures = df_in[input_col].map(Structure.from_dict).to_dict()

    # === Step 2. Relax Structures ===
    relax_results: dict[str, dict[str, Any]] = {}
    for material_id in tqdm(structures, desc=input_path.stem):
        structure, energy = relaxer.relax(
            structures[material_id],
            optimizer=optimizer,
            cell_filter=cell_filter,
            force_max=force_max,
            max_steps=max_steps,
        )
        relax_results[material_id] = {
            f"{model_name}_structure": structure,
            f"{model_name}_energy": energy,
        }

    # === Step 3. Save Shard Output ===
    df_out = pd.DataFrame(relax_results).T
    df_out.index.name = str(Key.mat_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.reset_index().to_json(
        output_path,
        default_handler=as_dict_handler,
        orient="records",
        lines=True,
    )


def build_parser() -> argparse.ArgumentParser:
    """
    Build the command line parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True, help="Directory containing WBM shard JSON files.")
    parser.add_argument("--checkpoint-path", type=Path, required=True, help="Path to the DPA4 DeePMD checkpoint.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for relaxed shard outputs.")
    parser.add_argument("--model-name", default="dpa4", help="Column prefix used in saved outputs.")
    parser.add_argument("--optimizer", choices=["FIRE", "LBFGS", "BFGS"], default="FIRE", help="ASE optimizer.")
    parser.add_argument("--cell-filter", choices=["frechet", "unit"], default="frechet", help="ASE cell filter.")
    parser.add_argument("--force-max", type=float, default=0.02, help="Force convergence threshold in eV/A.")
    parser.add_argument("--max-steps", type=int, default=500, help="Maximum optimization steps.")
    parser.add_argument("--num-jobs", type=int, default=1, help="Number of interleaved shard jobs.")
    parser.add_argument("--job-index", type=int, default=0, help="Index of this interleaved shard job.")
    parser.add_argument("--debug", action="store_true", help="Only run the first selected shard.")
    return parser


def run_relax(args: argparse.Namespace) -> None:
    """
    Run DPA4 WBM relaxations.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments.
    """
    if not args.checkpoint_path.is_file():
        raise FileNotFoundError(f"Checkpoint not found: {args.checkpoint_path}")
    if args.num_jobs < 1:
        raise ValueError("num_jobs must be positive")
    if not 0 <= args.job_index < args.num_jobs:
        raise ValueError("job_index must satisfy 0 <= job_index < num_jobs")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_params = {
        "data_dir": str(args.data_dir),
        "checkpoint_path": str(args.checkpoint_path),
        "versions": {dep: version(dep) for dep in ("deepmd-kit", "numpy", "torch")},
        str(Key.task_type): "IS2RE",
        "optimizer": args.optimizer,
        "cell_filter": args.cell_filter,
        "force_max": args.force_max,
        "max_steps": args.max_steps,
        "model_name": args.model_name,
        "num_jobs": args.num_jobs,
        "job_index": args.job_index,
    }
    with open(args.output_dir / f"run_params_{args.job_index}.json", mode="w") as file:
        json.dump(run_params, file, indent=4)

    shard_paths = collect_input_shards(args.data_dir)[args.job_index :: args.num_jobs]
    if args.debug:
        shard_paths = shard_paths[:1]

    relaxer = Relaxer(args.checkpoint_path)
    for input_path in shard_paths:
        output_path = args.output_dir / f"results_{input_path.stem}.json.gz"
        relax_shard(
            input_path=input_path,
            output_path=output_path,
            relaxer=relaxer,
            model_name=args.model_name,
            optimizer=args.optimizer,
            cell_filter=args.cell_filter,
            force_max=args.force_max,
            max_steps=args.max_steps,
        )


def main() -> None:
    """Run the DPA4 discovery CLI."""
    run_relax(build_parser().parse_args())


if __name__ == "__main__":
    main()
