import gzip
import traceback
import warnings
from copy import deepcopy
from io import BytesIO
from pathlib import Path
from typing import Any

import ase.io
import numpy as np
import pandas as pd
import requests  # type: ignore
import torch
import typer
import wandb
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from ase.spacegroup.symmetrize import check_symmetry
from orb_models.forcefield.calculator import ORBCalculator
from orb_models.forcefield.pretrained import ORB_PRETRAINED_MODELS
from tqdm import tqdm

from . import metrics
from . import thermal_conductivity as ltc
from .keys import Key

app = typer.Typer()


def download_results_json(data_root: Path) -> Path:
    """Download kappa results from figshare.

    url: https://figshare.com/ndownloader/files/52014656
    path: kappas-phononDB-PBE-noNAC.json
    """
    data_dir = data_root / "thermal_conductivity"
    data_dir.mkdir(exist_ok=True, parents=True)
    json_path = data_dir / "kappas-phononDB-PBE-noNAC.json"

    # Download and decompress in one step
    url = "https://figshare.com/ndownloader/files/52014656"
    response = requests.get(url)
    response.raise_for_status()

    # Decompress response content to json file
    with gzip.open(BytesIO(response.content)) as gz_file:
        with open(json_path, "wb") as json_file:
            json_file.write(gz_file.read())

    print(f"Downloaded and decompressed to {json_path}")
    return json_path


def download_phonondb_structures(data_root: Path) -> Path:
    """Download phonondb structures from figshare.

    url: https://figshare.com/ndownloader/files/51680888
    path: phonons/2024-11-09-phononDB-PBE-103-structures.extxyz
    """
    data_dir = data_root / "thermal_conductivity"
    data_dir.mkdir(exist_ok=True, parents=True)
    extxyz_path = data_dir / "2024-11-09-phononDB-PBE-103-structures.extxyz"

    url = "https://figshare.com/ndownloader/files/51680888"
    response = requests.get(url)
    response.raise_for_status()

    with open(extxyz_path, "wb") as file:
        file.write(response.content)

    print(f"Downloaded and saved to {extxyz_path}")
    return extxyz_path


def compute_and_log_final_results(
    df_kappa: pd.DataFrame,
    run: wandb.sdk.wandb_run.Run,
) -> None:
    """Compute final metrics and log results to Weights & Biases.

    Args:
        df_kappa: DataFrame containing thermal conductivity predictions
        run: Active W&B run for logging results
        job_name: Name of the job/experiment
        out_dir: Directory containing output files
    """
    json_path = download_results_json(Path("./data"))
    df_dft = pd.read_json(json_path).set_index("mp_id")

    df_ml_metrics = metrics.calc_kappa_metrics_from_dfs(df_kappa, df_dft)

    # log summary results
    kappa_sre = df_ml_metrics[Key.sre].mean()
    kappa_srme = df_ml_metrics[Key.srme].mean()
    run.summary.update({"srme": kappa_srme, "sre": kappa_sre})


def check_imaginary_freqs(frequencies: np.ndarray, threshold: float = -0.01) -> bool:
    """Check if frequencies are imaginary.

    Args:
        frequencies (np.ndarray): Frequencies to check
        threshold (float): Threshold for imaginary frequencies. Defaults to -0.01.
    Returns:
        bool: True if imaginary frequencies are found.
    """
    # Return True if all frequencies are NaN, indicating invalid or missing data
    if np.all(pd.isna(frequencies)):
        return True

    # Check for imaginary frequencies in non-acoustic modes at gamma point (q=0)
    # Indices 3+ correspond to optical modes which should never be negative
    if np.any(frequencies[0, 3:] < 0):
        return True

    # Check acoustic modes at gamma point against threshold. First 3 modes at q=0
    # are acoustic and may be slightly negative due to numerical noise
    if np.any(frequencies[0, :3] < threshold):
        return True

    # Check for imaginary frequencies at any q-point except gamma
    # All frequencies should be positive away from gamma point
    return bool(np.any(frequencies[1:] < 0))


def calc_kappa_for_structure(
    *,
    atoms: Atoms,
    calc: Calculator,
    displacement_distance: float,
    temperatures: list[float],
    ase_optimizer: str,
    ase_filter: str,
    max_steps: int,
    force_max: float,
    symprec: float,
    enforce_relax_symm: bool,
    ignore_broken_symm: bool,
    ignore_imaginary_freqs: bool,
    save_forces: bool,
    out_dir: str,
    task_id: int,
    is_plusminus: bool | str = True,
) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    """Predict ML kappa for single structure.

    Args:
        atoms (Atoms): Input structure
        calc (Calculator): ASE calculator to use
        displacement_distance (float): Displacement distance for phonon calculations
        temperatures (list[float]): Which temperatures to calculate kappa at in Kelvin
        ase_optimizer (str): ASE optimizer to use
        ase_filter (str): ASE filter to use
        max_steps (int): Maximum number of optimization steps
        force_max (float): Maximum force tolerance
        symprec (float): Symmetry precision
        enforce_relax_symm (bool): Whether to enforce symmetry during relaxation
        ignore_broken_symm (bool): If true, continue the calculation even though
            symmetry is broken.
        ignore_imaginary_freqs (bool): If true, continue the calculation even though
            imaginary frequencies are present.
        save_forces (bool): Whether to save force sets
        out_dir (str): Output directory
        task_id (int): Task ID for logging

    Returns:
        tuple[str, dict[str, Any], dict[str, Any] | None]:
            material ID, results dict, force results dict
    """
    warnings.filterwarnings("ignore", category=FutureWarning, module="torch")
    # Create a deep copy of the atoms
    atoms = atoms.copy()
    atoms.arrays = {k: v.copy() for k, v in atoms.arrays.items()}

    mat_id = atoms.info[Key.mat_id]
    init_info = deepcopy(atoms.info)
    mat_name = atoms.info["name"]
    info_dict = {
        "name": mat_name,
        "errors": [],
        "error_traceback": [],
    }

    ##################################
    #  Symmetry-preserving Relaxation
    ##################################
    filter_cls: type[ExpCellFilter | FrechetCellFilter] = {
        "frechet": FrechetCellFilter,
        "exp": ExpCellFilter,
    }[ase_filter]

    optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

    # Initialize variables that might be needed in error handling
    relax_dict = {"max_stress": None, "reached_max_steps": False}
    force_results = None
    try:
        atoms.calc = calc
        if enforce_relax_symm:
            atoms.set_constraint(FixSymmetry(atoms))
            filtered_atoms = filter_cls(atoms, mask=[True] * 3 + [False] * 3)
        else:
            filtered_atoms = filter_cls(atoms)

        optimizer = optim_cls(
            filtered_atoms,
            logfile=f"{out_dir}/relax_{task_id}.log",  # type: ignore
        )

        pre_symmetry_group = check_symmetry(atoms, symprec).number
        optimizer.run(fmax=force_max, steps=max_steps)
        post_symmetry_group = check_symmetry(atoms, symprec).number

        reached_max_steps = optimizer.step == max_steps
        if reached_max_steps:
            print(f"Material {mat_id=} reached {max_steps=} during relaxation")

        # maximum residual stress component in for xx,yy,zz and xy,yz,xz
        # components separately result is a array of 2 elements
        max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

        atoms.calc = None
        atoms.constraints = None
        atoms.info = init_info | atoms.info

        relax_dict = {
            "max_stress": max_stress,
            "reached_max_steps": reached_max_steps,
            "broken_symmetry": pre_symmetry_group != post_symmetry_group,
        }
        if not ignore_broken_symm and pre_symmetry_group != post_symmetry_group:
            raise ValueError(
                f"Symmetry group changed from {pre_symmetry_group} to {post_symmetry_group}"
            )

    except Exception as exc:
        warnings.warn(f"Failed to relax {mat_name=}, {mat_id=}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"RelaxError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict, None

    ####################
    #  Force Constants
    ####################
    try:
        ph3 = ltc.init_phono3py(
            atoms,
            fc2_supercell=atoms.info["fc2_supercell"],
            fc3_supercell=atoms.info["fc3_supercell"],
            q_point_mesh=atoms.info["q_mesh"],
            displacement_distance=displacement_distance,
            symprec=symprec,
            is_plusminus=is_plusminus,
        )

        ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
            ph3, calculator=calc, pbar_kwargs={"disable": True}
        )

        has_imaginary_freqs = check_imaginary_freqs(freqs)
        freqs_dict = {Key.has_imag_ph_modes: has_imaginary_freqs, Key.ph_freqs: freqs}

        continue_computing_conductivity = (
            not has_imaginary_freqs or ignore_imaginary_freqs
        )
        if continue_computing_conductivity:
            fc3_set = ltc.calculate_fc3_set(
                ph3, calculator=calc, pbar_kwargs={"position": task_id}
            )
            ph3.produce_fc3(symmetrize_fc3r=True)
        else:
            fc3_set = []  # type: ignore

        force_results = (
            {"fc2_set": fc2_set, "fc3_set": fc3_set} if save_forces else None
        )

        if not continue_computing_conductivity:
            warnings.warn(
                f"Skipping conductivity calculation for {mat_id} due to imaginary frequencies"
            )
            return mat_id, info_dict | relax_dict | freqs_dict, force_results

    except Exception as exc:
        warnings.warn(f"Failed to calculate force sets {mat_id}: {exc!r}", stacklevel=2)
        traceback.print_exc()
        info_dict["errors"].append(f"ForceConstantError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict, force_results

    ############################
    #  Conductivity Calculation
    ############################
    try:
        ph3, kappa_dict, _cond = ltc.calculate_conductivity(
            ph3, temperatures=temperatures
        )
        return mat_id, info_dict | relax_dict | freqs_dict | kappa_dict, force_results

    except Exception as exc:
        warnings.warn(
            f"Failed to calculate conductivity {mat_id}: {exc!r}", stacklevel=2
        )
        traceback.print_exc()
        info_dict["errors"].append(f"ConductivityError: {exc!r}")
        info_dict["error_traceback"].append(traceback.format_exc())
        return mat_id, info_dict | relax_dict | freqs_dict, force_results


@app.command()
def run_thermal_conductivity(
    model_name_or_uri: str = "orb-v3-conservative-inf-mpa",  # or orb-v3-conservative-inf-omat
    limit_structures: int | None = None,
    displacement_distance: float = 0.03,
    ase_optimizer: str = "FIRE",
    ase_filter: str = "frechet",
    max_steps: int = 300,
    force_max: float = 1e-4,  # Run until the forces are smaller than this in eV/A
    symprec: float = 1e-5,  # symmetry precision for enforcing relaxation and conductivity calcs
    enforce_relax_symm: bool = True,  # Enforce symmetry during relaxation if broken
    ignore_broken_symm: bool = False,  # Continue calculation even if symmetry group changed during relaxation
    ignore_imaginary_freqs: bool = False,  # Continue calculation even if imaginary freqs are present in fc2
    save_forces: bool = True,  # Save force sets to file
    temperatures: list[float] = [300],
    precision: str = "float64",
    deterministic: bool = False,
    max_num_neighbors: int | None = 120,
    wandb_name: str | None = None,
    is_plusminus: bool = True,  # False corresponds to "auto" (None would be better, but typer won't allow)
) -> None:
    """Run thermal conductivity prediction benchmark, as used by matbench-discovery.

    This script:
    1. Downloads structures from PhononDB
    2. Relaxes each structure using the specified model
    3. Calculates force constants and thermal conductivity
    4. Saves intermediate quantities (like force sets)
    5. Logs SRME and SRE metrics to wandb
    """
    model_name = model_name_or_uri.split("/")[-1].replace(":", "-")
    project_name = (
        "thermal-conductivity"
        if limit_structures is None
        else "thermal-conductivity-test"
    )
    run_name = f"{wandb_name or model_name}-disp{displacement_distance}-{precision}"
    if deterministic:
        run_name += "-deterministic"
    if not is_plusminus:
        run_name += "-auto"
    is_plusminus = is_plusminus if is_plusminus else "auto"  # type: ignore

    run_params = {
        "model_name_or_uri": model_name_or_uri,
        "limit_structures": limit_structures,
        "displacement_distance": displacement_distance,
        "is_plusminus": is_plusminus,
        "max_num_neighbors": max_num_neighbors,
        "precision": precision,
    }

    run = wandb.init(project=project_name, name=run_name, config=run_params)

    out_dir = Path(run_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Results will be saved to {out_dir}")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using {device=}")

    model = ORB_PRETRAINED_MODELS[model_name_or_uri]()
    model.to(device)
    calc = ORBCalculator(
        model,
        max_num_neighbors=max_num_neighbors,
        device=device,
    )
    if deterministic:
        torch.use_deterministic_algorithms(True)

    structure_path = download_phonondb_structures(Path("./data"))
    atoms_list = ase.io.read(structure_path, index=":")
    atoms_list = atoms_list[:limit_structures] if limit_structures else atoms_list

    kappa_results: dict[str, dict[str, Any]] = {}
    force_results: dict[str, dict[str, Any]] = {}
    for idx, atoms in tqdm(
        enumerate(atoms_list), desc=f"Predicting kappa with {model_name_or_uri}"
    ):
        mat_id, result_dict, force_dict = calc_kappa_for_structure(
            atoms=atoms,
            calc=calc,
            displacement_distance=displacement_distance,
            temperatures=temperatures,
            ase_optimizer=ase_optimizer,
            ase_filter=ase_filter,
            max_steps=max_steps,
            force_max=force_max,
            symprec=symprec,
            enforce_relax_symm=enforce_relax_symm,
            ignore_broken_symm=ignore_broken_symm,
            ignore_imaginary_freqs=ignore_imaginary_freqs,
            save_forces=save_forces,
            out_dir=str(out_dir),
            task_id=idx,
            is_plusminus=is_plusminus,
        )
        kappa_results[mat_id] = result_dict
        if force_dict is not None:
            force_results[mat_id] = force_dict

    df_kappa = pd.DataFrame(kappa_results).T
    df_kappa.index.name = Key.mat_id
    df_kappa.reset_index().to_json(f"{out_dir}/kappa.json.gz")

    if save_forces:
        df_force = pd.DataFrame(force_results).T
        df_force = pd.concat([df_kappa, df_force], axis=1)
        df_force.index.name = Key.mat_id
        df_force.reset_index().to_json(f"{out_dir}/force-sets.json.gz")

    # Calculate final SRME and SRE metrics
    compute_and_log_final_results(df_kappa, run)


if __name__ == "__main__":
    app()
