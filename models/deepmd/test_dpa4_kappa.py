"""Run DPA4-series thermal-conductivity (kappa) calculations for Matbench Discovery."""

from __future__ import annotations

import argparse
import json
import traceback
import warnings
from copy import deepcopy
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
import torch
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from deepmd.calculator import DP
from moyopy import MoyoDataset
from moyopy.interface import MoyoAdapter
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.enums import DataFiles, MbdKey
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons import thermal_conductivity as ltc

if TYPE_CHECKING:
    from collections.abc import Callable

    from ase import Atoms
    from ase.optimize.optimize import Optimizer

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")

OPTIM_CLS: dict[str, Callable[..., Optimizer]] = {"FIRE": FIRE, "LBFGS": LBFGS}

# Conductivity tensor columns that phono3py reports in 6-component Voigt notation
# and that must be expanded to full 3x3 tensors for the Matbench evaluator.
KAPPA_TENSOR_KEYS = (
    MbdKey.kappa_tot_rta,
    MbdKey.kappa_p_rta,
    MbdKey.kappa_c,
    MbdKey.mode_kappa_tot_rta,
)


def voigt_6_to_full_3x3(tensor: np.ndarray) -> np.ndarray:
    """Expand a trailing Voigt-6 axis (xx, yy, zz, yz, xz, xy) to a full 3x3 matrix.

    Depending on the phono3py version, conductivity tensors are reported either as
    full (..., 3, 3) matrices or in 6-component Voigt notation (..., 6). The Matbench
    Discovery kappa evaluator (``calculate_kappa_avg``) requires a trailing (3, 3) or
    (3,) axis, so Voigt-packed tensors must be expanded before saving. Arrays whose
    trailing axis is not length 6 are returned unchanged (already full / not a tensor).
    """
    arr = np.asarray(tensor, dtype=float)
    if arr.ndim == 0 or arr.shape[-1] != 6:
        return tensor
    xx, yy, zz, yz, xz, xy = (arr[..., idx] for idx in range(6))
    return np.stack(
        [
            np.stack([xx, xy, xz], axis=-1),
            np.stack([xy, yy, yz], axis=-1),
            np.stack([xz, yz, zz], axis=-1),
        ],
        axis=-2,
    )


class KappaRunner:
    """DPA4 kappa benchmark runner using the Matbench phonon workflow."""

    def __init__(
        self,
        checkpoint_path: Path,
        output_dir: Path,
        structures_path: Path | None = None,
        optimizer: Literal["FIRE", "LBFGS"] = "FIRE",
        max_steps: int = 300,
        force_max: float = 1e-4,
        symprec: float = 1e-5,
        displacement_distance: float = 0.03,
        limit: int = 0,
        *,
        save_forces: bool = True,
        enforce_relax_symm: bool = True,
        conductivity_broken_symm: bool = False,
        prog_bar: bool = True,
    ) -> None:
        """Initialize benchmark settings.

        Parameters
        ----------
        checkpoint_path : Path
            Path to the DPA4 DeePMD checkpoint.
        output_dir : Path
            Directory for benchmark outputs.
        structures_path : Path | None, default=None
            Optional phononDB structure file. Defaults to Matbench Discovery data.
        optimizer : {"FIRE", "LBFGS"}, default="FIRE"
            ASE optimizer used for relaxation.
        max_steps : int, default=300
            Maximum relaxation steps.
        force_max : float, default=1e-4
            Force convergence threshold in eV/A.
        symprec : float, default=1e-5
            Symmetry precision for relaxation and conductivity.
        displacement_distance : float, default=0.03
            Phono3py displacement distance in Angstrom.
        limit : int, default=0
            Number of structures to run. Use 0 for all structures.
        save_forces : bool, default=True
            Whether to save FC2/FC3 force sets.
        enforce_relax_symm : bool, default=True
            Whether to enforce symmetry during relaxation.
        conductivity_broken_symm : bool, default=False
            Whether to calculate conductivity after symmetry breaking.
        prog_bar : bool, default=True
            Whether to show progress bars.
        """
        self.checkpoint_path = checkpoint_path
        self.output_dir = output_dir
        self.structures_path = structures_path
        self.optimizer = optimizer
        self.max_steps = max_steps
        self.force_max = force_max
        self.symprec = symprec
        self.displacement_distance = displacement_distance
        self.limit = limit
        self.save_forces = save_forces
        self.enforce_relax_symm = enforce_relax_symm
        self.conductivity_broken_symm = conductivity_broken_symm
        self.prog_bar = prog_bar

    def run(self) -> Path:
        """Run the kappa benchmark and save outputs.

        Returns:
        -------
        Path
            Path to the saved conductivity JSON file.
        """
        # === Step 1. Initialize Model and Inputs ===
        self.output_dir.mkdir(parents=True, exist_ok=True)
        calculator = DP(str(self.checkpoint_path))
        optim_cls = OPTIM_CLS[self.optimizer]
        structures_path = self.structures_path or Path(
            DataFiles.phonondb_pbe_103_structures.path
        )
        atoms_list = read(structures_path, format="extxyz", index=":")
        if self.limit > 0:
            atoms_list = atoms_list[: self.limit]

        timestamp = f"{datetime.now().astimezone():%Y-%m-%d %H:%M:%S}"
        run_params = {
            "timestamp": timestamp,
            "checkpoint_path": str(self.checkpoint_path),
            "structures_path": str(structures_path),
            "device": "cuda" if torch.cuda.is_available() else "cpu",
            "versions": {
                dep: version(dep)
                for dep in ("deepmd-kit", "numpy", "torch", "matbench_discovery")
            },
            "ase_optimizer": self.optimizer,
            "cell_filter": "FrechetCellFilter",
            "max_steps": self.max_steps,
            "force_max": self.force_max,
            "symprec": self.symprec,
            "enforce_relax_symm": self.enforce_relax_symm,
            "conductivity_broken_symm": self.conductivity_broken_symm,
            "temperatures": [300],
            "displacement_distance": self.displacement_distance,
            "n_structures": len(atoms_list),
        }
        with open(self.output_dir / "run_params.json", mode="w") as file:
            json.dump(run_params, file, indent=4)

        # === Step 2. Relax and Calculate Conductivity ===
        force_results: dict[str, dict[str, Any]] = {}
        kappa_results: dict[str, dict[str, Any]] = {}
        progress = tqdm(
            atoms_list,
            desc="Conductivity calculation",
            disable=not self.prog_bar,
        )

        for atoms in progress:
            self._run_one_structure(
                atoms=atoms,
                calculator=calculator,
                optim_cls=optim_cls,
                force_results=force_results,
                kappa_results=kappa_results,
                progress=progress,
            )

        # === Step 3. Save Outputs ===
        out_path = self.output_dir / "conductivity.json.gz"
        df_kappa = pd.DataFrame(kappa_results).T
        df_kappa.index.name = str(Key.mat_id)
        df_kappa.reset_index().to_json(out_path, compression="gzip")
        print(f"Saved kappa results to {out_path}")

        if self.save_forces:
            force_out_path = self.output_dir / "force_sets.json.gz"
            df_force = pd.DataFrame(force_results).T
            df_force.index.name = str(Key.mat_id)
            df_force.reset_index().to_json(force_out_path, compression="gzip")
            print(f"Saved force results to {force_out_path}")

        return out_path

    def _run_one_structure(
        self,
        atoms: Atoms,
        calculator: DP,
        optim_cls: Callable[..., Optimizer],
        force_results: dict[str, dict[str, Any]],
        kappa_results: dict[str, dict[str, Any]],
        progress: tqdm,
    ) -> None:
        """Run relaxation and conductivity for one structure.

        Parameters
        ----------
        atoms : Atoms
            Input phononDB structure.
        calculator : DP
            DeePMD calculator.
        optim_cls : Callable[..., Optimizer]
            ASE optimizer class.
        force_results : dict[str, dict[str, Any]]
            Accumulator for force sets.
        kappa_results : dict[str, dict[str, Any]]
            Accumulator for conductivity results.
        progress : tqdm
            Progress bar to update.
        """
        mat_id = atoms.info.get(str(Key.mat_id), atoms.info.get("mp_id"))
        if mat_id is None:
            raise KeyError(
                f"No material ID found in atoms.info keys: {sorted(atoms.info)}"
            )
        init_info = deepcopy(atoms.info)
        formula = atoms.get_chemical_formula()
        spg_num = MoyoDataset(
            MoyoAdapter.from_atoms(atoms), symprec=self.symprec
        ).number
        info_dict: dict[str, Any] = {
            str(Key.formula): formula,
            str(Key.spg_num): spg_num,
        }
        err_dict: dict[str, list[str]] = {"errors": [], "error_traceback": []}
        progress.set_postfix_str(mat_id, refresh=True)

        relax_dict: dict[str, Any] = {
            "max_stress": None,
            "reached_max_steps": False,
            "broken_symmetry": False,
        }

        try:
            relax_dict = self._relax_atoms(
                atoms=atoms,
                calculator=calculator,
                optim_cls=optim_cls,
                init_info=init_info,
                spg_num=spg_num,
                mat_id=mat_id,
            )
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            warnings.warn(
                f"Failed to relax {formula=}, {mat_id=}: {exc!r}", stacklevel=2
            )
            traceback.print_exc()
            err_dict["errors"].append(f"RelaxError: {exc!r}")
            err_dict["error_traceback"].append(traceback.format_exc())
            kappa_results[mat_id] = info_dict | relax_dict | err_dict
            return

        try:
            q_point_mesh = atoms.info.get("q_point_mesh", atoms.info.get("q_mesh"))
            if q_point_mesh is None:
                raise ValueError(
                    f"missing q_point_mesh/q_mesh keys in atoms.info for {mat_id}"
                )
            ph3 = ltc.init_phono3py(
                atoms,
                fc2_supercell=atoms.info["fc2_supercell"],
                fc3_supercell=atoms.info["fc3_supercell"],
                q_point_mesh=q_point_mesh,
                displacement_distance=self.displacement_distance,
                symprec=self.symprec,
            )
            ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
                ph3,
                calculator=calculator,
                pbar_kwargs={"leave": False, "disable": not self.prog_bar},
            )
            freqs_dict = {
                Key.has_imag_ph_modes: check_imaginary_freqs(freqs),
                Key.ph_freqs: freqs,
            }

            run_ltc = not freqs_dict[Key.has_imag_ph_modes] and (
                not relax_dict["broken_symmetry"] or self.conductivity_broken_symm
            )
            if run_ltc:
                print(f"Calculating FC3 for {mat_id}")
                fc3_set = ltc.calculate_fc3_set(
                    ph3,
                    calculator=calculator,
                    pbar_kwargs={"leave": False, "disable": not self.prog_bar},
                )
                ph3.produce_fc3(symmetrize_fc3r=True)
            else:
                fc3_set = []

            if self.save_forces:
                force_results[mat_id] = {"fc2_set": fc2_set, "fc3_set": fc3_set}
            if not run_ltc:
                kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | err_dict
                return
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            warnings.warn(
                f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2
            )
            traceback.print_exc()
            err_dict["errors"].append(f"ForceConstantError: {exc!r}")
            err_dict["error_traceback"].append(traceback.format_exc())
            kappa_results[mat_id] = info_dict | relax_dict | err_dict
            return

        try:
            ph3, kappa_dict, _cond = ltc.calculate_conductivity(ph3, temperatures=[300])
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            warnings.warn(
                f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
            )
            traceback.print_exc()
            err_dict["errors"].append(f"ConductivityError: {exc!r}")
            err_dict["error_traceback"].append(traceback.format_exc())
            kappa_results[mat_id] = info_dict | relax_dict | freqs_dict | err_dict
            return

        # phono3py may pack conductivity tensors in 6-component Voigt notation;
        # expand them to full 3x3 so the saved file matches the layout the Matbench
        # Discovery kappa evaluator expects (trailing (3, 3) axis).
        for kappa_key in KAPPA_TENSOR_KEYS:
            if kappa_key in kappa_dict:
                kappa_dict[kappa_key] = voigt_6_to_full_3x3(kappa_dict[kappa_key])

        kappa_results[mat_id] = (
            info_dict | relax_dict | freqs_dict | kappa_dict | err_dict
        )

    def _relax_atoms(
        self,
        atoms: Atoms,
        calculator: DP,
        optim_cls: Callable[..., Optimizer],
        init_info: dict[str, Any],
        spg_num: int,
        mat_id: str,
    ) -> dict[str, Any]:
        """Relax one structure before force-constant calculation.

        Parameters
        ----------
        atoms : Atoms
            Input phononDB structure.
        calculator : DP
            DeePMD calculator.
        optim_cls : Callable[..., Optimizer]
            ASE optimizer class.
        init_info : dict[str, Any]
            Original atoms info dictionary.
        spg_num : int
            Initial space group number.
        mat_id : str
            Material identifier used for log files.

        Returns:
        -------
        dict[str, Any]
            Relaxation metadata.
        """
        atoms.calc = calculator
        if self.max_steps <= 0:
            return {
                "max_stress": None,
                "reached_max_steps": False,
                "broken_symmetry": False,
            }

        if self.enforce_relax_symm:
            atoms.set_constraint(FixSymmetry(atoms))
        filtered_atoms = FrechetCellFilter(atoms, mask=[True] * 3 + [False] * 3)
        optimizer = optim_cls(
            filtered_atoms, logfile=self.output_dir / f"relax_{mat_id}.log"
        )
        optimizer.run(fmax=self.force_max, steps=self.max_steps)

        reached_max_steps = optimizer.nsteps >= self.max_steps
        if reached_max_steps:
            print(f"{mat_id=} reached {self.max_steps=}")

        max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)
        atoms.calc = None
        atoms.constraints = None
        atoms.info = init_info | atoms.info

        relaxed_spg = MoyoDataset(
            MoyoAdapter.from_atoms(atoms), symprec=self.symprec
        ).number
        return {
            "max_stress": max_stress,
            "reached_max_steps": reached_max_steps,
            "relaxed_space_group_number": relaxed_spg,
            "broken_symmetry": spg_num != relaxed_spg,
        }


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
    -------
    argparse.ArgumentParser
        Configured parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--checkpoint-path",
        type=Path,
        required=True,
        help="Path to the DPA4 DeePMD checkpoint.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for kappa files.",
    )
    parser.add_argument(
        "--structures-path",
        type=Path,
        default=None,
        help="Optional phononDB extxyz path.",
    )
    parser.add_argument(
        "--optimizer", choices=["FIRE", "LBFGS"], default="FIRE", help="ASE optimizer."
    )
    parser.add_argument(
        "--max-steps", type=int, default=300, help="Maximum relaxation steps."
    )
    parser.add_argument(
        "--force-max",
        type=float,
        default=1e-4,
        help="Force convergence threshold in eV/A.",
    )
    parser.add_argument(
        "--symprec", type=float, default=1e-5, help="Symmetry precision."
    )
    parser.add_argument(
        "--displacement-distance",
        type=float,
        default=0.03,
        help="Phono3py displacement distance in Angstrom.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Number of structures to run. Use 0 for all.",
    )
    parser.add_argument(
        "--no-save-forces", action="store_true", help="Skip writing force_sets.json.gz."
    )
    return parser


def kappa_run(args: argparse.Namespace) -> None:
    """Run the DPA4 kappa benchmark.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments.
    """
    if not args.checkpoint_path.is_file():
        raise FileNotFoundError(f"Checkpoint not found: {args.checkpoint_path}")
    if args.structures_path is not None and not args.structures_path.is_file():
        raise FileNotFoundError(f"Structure file not found: {args.structures_path}")
    if args.limit < 0:
        raise ValueError("limit must be non-negative")

    runner = KappaRunner(
        checkpoint_path=args.checkpoint_path,
        output_dir=args.output_dir,
        structures_path=args.structures_path,
        optimizer=args.optimizer,
        max_steps=args.max_steps,
        force_max=args.force_max,
        symprec=args.symprec,
        displacement_distance=args.displacement_distance,
        limit=args.limit,
        save_forces=not args.no_save_forces,
    )
    runner.run()


def main() -> None:
    """Run the DPA4 kappa CLI."""
    kappa_run(build_parser().parse_args())


if __name__ == "__main__":
    main()
