# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy>=1.23",
#   "pandas>=2.0",
#   "torch>=2.0",
#   "ase>=3.22",
#   "pymatgen>=2024.0",
#   "pymatviz>=0.15",
#   "tqdm>=4.65",
#   "matbench-discovery",
#   "bam-torch",
# ]
# ///
"""WBM IS2RE-SR evaluation for BAM-MP-core (matbench-discovery leaderboard).

Relaxes every WBM initial structure with ASE's FIRE optimizer wrapped in
``FrechetCellFilter`` (cell + positions) using the BAM-torch ASE-compatible
``RACECalculator``. Produces the discovery CSV and geo-opt JSONL files
declared in ``models/BAM/bam-mp-core.yml``.

Multi-GPU parallelisation is done by setting two environment variables and
launching N copies of this script, one per visible GPU:

    CUDA_VISIBLE_DEVICES=<id>  RACE_N_SPLITS=<N>  RACE_SPLIT_ID=<k>  python ...

Outputs (relative to this file):

    ./results/{today}-wbm-IS2RE-FIRE/
        results-000.json.gz      # per-split relaxed structures + energies
        checkpoint-000.json.gz   # resumable checkpoint (deleted on completion)
        relax_log-000.tsv        # per-structure step/fmax/time log
        traj-000/{mat_id}.traj   # full ASE trajectories (optional, large)

After all splits complete, merge the per-split JSONL files, convert relaxed
energies to formation energies, and write the artifacts declared in
``bam-mp-core.yml``.
"""

from __future__ import annotations

import os
import shutil
import time
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Protocol
from urllib.request import urlopen

import numpy as np
import pandas as pd
import torch
from ase.calculators.calculator import all_changes
from ase.filters import FrechetCellFilter
from ase.io import Trajectory
from ase.optimize import FIRE
from bam_torch.tase.base_calculator import RACECalculator as BaseRACECalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import DataFiles

if TYPE_CHECKING:
    from ase import Atoms

__author__ = "Hyun Gyu Park"
__date__ = "2026-04-22"


# -----------------------------------------------------------------------------
# Configuration — override via environment variables
# -----------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent

# Model checkpoint: set BAM_MODEL_PKL to a local checkpoint path. If that file is
# absent, the script downloads BAM_MODEL_PKL_URL into the same cache path. Replace
# TODO_UPLOAD_DIRECT_URL with the Figshare direct-download URL after upload.
BAM_MODEL_PKL_URL = os.environ.get(
    "BAM_MODEL_PKL_URL",
    "https://ndownloader.figshare.com/files/65859867",
)
CKPT_PATH = Path(
    os.environ.get(
        "BAM_MODEL_PKL",
        Path.home() / ".cache" / "bam_torch" / "ckpt_race_best.pkl",
    )
).expanduser()


def ensure_checkpoint(path: Path) -> Path:
    if path.is_file():
        print(f"Using BAM-MP-core checkpoint: {path}")
        return path
    if not BAM_MODEL_PKL_URL or BAM_MODEL_PKL_URL.startswith("TODO"):
        raise FileNotFoundError(
            f"BAM checkpoint not found: {path}. Set BAM_MODEL_PKL to a local "
            "checkpoint path or set BAM_MODEL_PKL_URL to a direct download URL."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading BAM-MP-core checkpoint to {path}...")
    tmp_path = path.with_suffix(path.suffix + ".part")
    try:
        with (
            urlopen(BAM_MODEL_PKL_URL, timeout=60) as response,  # noqa: S310
            tmp_path.open("wb") as file,
        ):
            shutil.copyfileobj(response, file)
            expected_size = response.headers.get("Content-Length")
        size = tmp_path.stat().st_size
        if size == 0:
            raise OSError(f"Downloaded checkpoint is empty: {tmp_path}")
        if expected_size is not None and size != int(expected_size):
            raise OSError(
                f"Incomplete checkpoint download: expected {expected_size} bytes, "
                f"got {size}"
            )
        tmp_path.replace(path)
    finally:
        tmp_path.unlink(missing_ok=True)
    print("Download complete.")
    return path


# Inference-time ASE/FIRE hyperparameters (mirrored in bam-mp-core.yml `hyperparams`)
CUTOFF = 6.0
FMAX = 0.05  # eV/Å — FIRE convergence threshold
MAX_STEPS = 500  # FIRE step cap
HEAD_IDX = None  # single-head pretrain checkpoint
SAVE_EVERY = 500  # checkpoint every N successful structures

# Multi-GPU splitting (set by the launcher shell script)
N_SPLITS = int(os.getenv("RACE_N_SPLITS", "1"))
SPLIT_ID = int(os.getenv("RACE_SPLIT_ID", "0"))
MODEL_TAG = "bam-mp-core"


# -----------------------------------------------------------------------------
# Load BAM-torch RACECalculator for BAM-MP-core
# -----------------------------------------------------------------------------

# Bug workaround: bam_torch.utils.data.batch_graphs assumes globals["stress"]
# exists when periodic=True, but preprocess_graph(targets=False) on the
# inference path does not populate it. Patch before the calculator loads.
import bam_torch.utils.data as _bt_data  # noqa: E402

_orig_batch_graphs = _bt_data.batch_graphs


class _BatchGraph(Protocol):
    globals: dict[str, np.ndarray]


def _batch_graphs_safe(graphs: list[_BatchGraph], periodic: bool = True) -> object:  # noqa: FBT001, FBT002
    if periodic:
        for g in graphs:
            if "stress" not in g.globals:
                g.globals["stress"] = np.zeros((1, 6), dtype=np.float32)
    return _orig_batch_graphs(graphs, periodic=periodic)


_bt_data.batch_graphs = _batch_graphs_safe


class BAMCalculator(BaseRACECalculator):
    """RACECalculator wrapper exposing free_energy for ASE cell filters."""

    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "forces",
        "stress",
        "free_energy",
    ]

    def calculate(
        self,
        atoms: Atoms,
        properties: list[str] | tuple[str, ...] = ("energy",),
        system_changes: list[str] | None = None,
    ) -> None:  # type: ignore[override]
        if system_changes is None:
            system_changes = all_changes
        super().calculate(atoms, list(properties), system_changes)
        if "energy" in self.results:
            self.results["free_energy"] = self.results["energy"]


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------
ARTIFACT_DIR = SCRIPT_DIR / f"{MODEL_TAG}-1.0v"
DISCOVERY_ARTIFACT = ARTIFACT_DIR / "2026-05-09-wbm-IS2RE.csv.gz"
GEO_OPT_ARTIFACT = ARTIFACT_DIR / "2026-05-09-wbm-geo-opt-FIRE.jsonl.gz"
PRED_COL = f"e_form_per_atom_{MODEL_TAG}"
STRUCT_COL = f"{MODEL_TAG}_structure"
ENERGY_COL = f"{MODEL_TAG}_energy"


def save_checkpoint(results: dict, path: Path) -> None:
    if not results:
        return
    df = pd.DataFrame(results).T.add_prefix(f"{MODEL_TAG}_")
    df.index.name = Key.mat_id
    df.reset_index().to_json(
        path, default_handler=as_dict_handler, orient="records", lines=True
    )


def _result_paths(out_dir: Path, out_path: Path) -> list[Path] | None:
    if N_SPLITS <= 1:
        return [out_path] if out_path.exists() else None

    paths = [out_dir / f"results-{idx:>03}.json.gz" for idx in range(N_SPLITS)]
    missing = [path.name for path in paths if not path.exists()]
    if missing:
        print(
            "Skipping registered artifact export until all split outputs exist: "
            + ", ".join(missing)
        )
        return None
    return paths


def _load_relaxation_results(result_paths: list[Path]) -> pd.DataFrame:
    df_results = pd.concat(
        [pd.read_json(path, orient="records", lines=True) for path in result_paths],
        ignore_index=True,
    )
    mat_id_col = str(Key.mat_id)
    return df_results.drop_duplicates(subset=mat_id_col, keep="last").set_index(
        mat_id_col
    )


def _formation_energy(row: pd.Series) -> float:
    from pymatgen.core import Structure

    from matbench_discovery.energy import (
        calc_energy_from_e_refs,
        mp_elemental_ref_energies,
    )

    structure = row[STRUCT_COL]
    if isinstance(structure, dict):
        structure = Structure.from_dict(structure)
    return calc_energy_from_e_refs(
        structure,
        ref_energies=mp_elemental_ref_energies,
        total_energy=float(row[ENERGY_COL]),
    )


def write_registered_artifacts(out_dir: Path, out_path: Path) -> None:
    result_paths = _result_paths(out_dir, out_path)
    if result_paths is None:
        return

    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    df_results = _load_relaxation_results(result_paths)
    df_results.reset_index().to_json(
        GEO_OPT_ARTIFACT,
        default_handler=as_dict_handler,
        orient="records",
        lines=True,
    )

    df_preds = pd.DataFrame({PRED_COL: df_results.apply(_formation_energy, axis=1)})
    df_preds.index.name = Key.mat_id
    df_preds[PRED_COL].round(4).to_csv(DISCOVERY_ARTIFACT)
    print(
        f"Registered artifacts written:\n  {DISCOVERY_ARTIFACT}\n  {GEO_OPT_ARTIFACT}"
    )


# -----------------------------------------------------------------------------
# Benchmark entry point
# -----------------------------------------------------------------------------
def main() -> None:
    ckpt_path = ensure_checkpoint(CKPT_PATH)

    out_dir = SCRIPT_DIR / "results" / f"{today}-wbm-IS2RE-FIRE"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (
        f"results-{SPLIT_ID:>03}.json.gz" if N_SPLITS > 1 else "results.json.gz"
    )
    ckpt_save_path = out_dir / (
        f"checkpoint-{SPLIT_ID:>03}.json.gz" if N_SPLITS > 1 else "checkpoint.json.gz"
    )

    # -------------------------------------------------------------------------
    # Load WBM initial structures (auto-downloaded via DataFiles)
    # -------------------------------------------------------------------------
    wbm_zip_path = Path(DataFiles.wbm_initial_atoms.path)
    print(f"Loading WBM structures from {wbm_zip_path} ...")
    atoms_list = ase_atoms_from_zip(wbm_zip_path)
    print(f"Loaded {len(atoms_list):,} structures total")

    if N_SPLITS > 1:
        total = len(atoms_list)
        chunk = total // N_SPLITS
        start = SPLIT_ID * chunk
        end = total if SPLIT_ID == N_SPLITS - 1 else start + chunk
        atoms_list = atoms_list[start:end]
        print(f"Split {SPLIT_ID}/{N_SPLITS}: {len(atoms_list):,} of {total:,}")

    # -------------------------------------------------------------------------
    # Resume from checkpoint if present
    # -------------------------------------------------------------------------
    relax_results: dict[str, dict] = {}
    if ckpt_save_path.exists():
        try:
            df_ckpt = pd.read_json(ckpt_save_path, orient="records", lines=True)
            for _, row in df_ckpt.iterrows():
                mat_id = row[str(Key.mat_id)]
                relax_results[mat_id] = {
                    "structure": row[f"{MODEL_TAG}_structure"],
                    "energy": row[f"{MODEL_TAG}_energy"],
                }
            print(f"Resumed from checkpoint: {len(relax_results):,} structures done")
        except Exception as exc:  # noqa: BLE001
            print(f"Could not load checkpoint ({exc!r}); starting fresh.")
            relax_results = {}

    # -------------------------------------------------------------------------
    # Load checkpoint-backed calculator
    # -------------------------------------------------------------------------
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    race_calc = BAMCalculator(model=str(ckpt_path), device=device)

    # -------------------------------------------------------------------------
    # Relaxation loop
    # -------------------------------------------------------------------------
    traj_dir = out_dir / (f"traj-{SPLIT_ID:>03}" if N_SPLITS > 1 else "traj")
    traj_dir.mkdir(parents=True, exist_ok=True)
    relax_log_path = out_dir / (
        f"relax_log-{SPLIT_ID:>03}.tsv" if N_SPLITS > 1 else "relax_log.tsv"
    )
    log_header_needed = not relax_log_path.exists()
    n_success = len(relax_results)
    n_failed = 0
    n_attempted = 0
    consecutive_failures = 0
    max_consecutive_failures = 10
    max_initial_failures = 25
    relax_times: list[float] = []
    total_start = time.time()

    relax_log_f = relax_log_path.open("a", buffering=1)
    try:
        if log_header_needed:
            relax_log_f.write(
                "mat_id\tn_atoms\tn_steps\tfinal_fmax\tconverged\telapsed_s\n"
            )

        for atoms_orig in tqdm(atoms_list, desc="Relaxing"):
            mat_id = atoms_orig.info[Key.mat_id]
            if mat_id in relax_results:
                continue
            atoms = deepcopy(atoms_orig)
            n_attempted += 1
            try:
                t0 = time.time()
                atoms.calc = race_calc
                traj_path = traj_dir / f"{mat_id}.traj"
                traj_writer = Trajectory(str(traj_path), "w", atoms)
                try:
                    dyn = FIRE(FrechetCellFilter(atoms), logfile=None)
                    dyn.attach(traj_writer)
                    traj_writer.write()
                    dyn.run(fmax=FMAX, steps=MAX_STEPS)
                finally:
                    traj_writer.close()

                n_steps = dyn.get_number_of_steps()
                final_fmax = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
                dt = time.time() - t0
                relax_log_f.write(
                    f"{mat_id}\t{len(atoms)}\t{n_steps}\t{final_fmax:.4f}\t"
                    f"{n_steps < MAX_STEPS}\t{dt:.2f}\n"
                )

                relax_results[mat_id] = {
                    "structure": AseAtomsAdaptor.get_structure(atoms),
                    "energy": atoms.get_potential_energy(),
                }
                relax_times.append(dt)
                n_success += 1
                consecutive_failures = 0

                if n_success % 100 == 0:
                    avg_t = sum(relax_times[-500:]) / min(len(relax_times), 500)
                    elapsed = time.time() - total_start
                    remaining = avg_t * (len(atoms_list) - n_success - n_failed)
                    print(
                        f"  [{n_success}/{len(atoms_list)}] avg={avg_t:.2f}s/struct, "
                        f"elapsed={elapsed:.0f}s, remaining~{remaining:.0f}s"
                    )
                if n_success % SAVE_EVERY == 0:
                    save_checkpoint(relax_results, ckpt_save_path)
                    print(f"  Checkpoint saved: {n_success} structures")
            except Exception as exc:
                n_failed += 1
                consecutive_failures += 1
                if n_failed <= 3:
                    import traceback

                    traceback.print_exc()
                else:
                    print(f"Failed {mat_id}: {exc!r}")
                if consecutive_failures >= max_consecutive_failures:
                    raise RuntimeError(
                        "Aborting BAM-MP-core relaxation after "
                        f"{consecutive_failures} consecutive failures."
                    ) from exc
                if n_success == 0 and n_failed >= max_initial_failures:
                    raise RuntimeError(
                        "Aborting BAM-MP-core relaxation after "
                        f"{n_failed}/{n_attempted} initial failures and no successes."
                    ) from exc
    finally:
        relax_log_f.close()

    # -------------------------------------------------------------------------
    # Finalize split and registered artifacts
    # -------------------------------------------------------------------------
    total_time = time.time() - total_start
    if not relax_results:
        raise RuntimeError(
            "No successful BAM-MP-core relaxations; refusing to write artifacts."
        )
    save_checkpoint(relax_results, out_path)
    if ckpt_save_path.exists():
        ckpt_save_path.unlink()
        print("Checkpoint removed (run complete)")
    write_registered_artifacts(out_dir, out_path)

    relax_times_arr = np.array(relax_times) if relax_times else np.array([0.0])
    bar = "=" * 60
    print(
        f"\n{bar}\nSplit {SPLIT_ID} Complete\n{bar}\n"
        f"  Structures: {n_success} success, {n_failed} failed\n"
        f"  Total time: {total_time:.1f}s ({total_time / 3600:.2f}h)\n"
        f"  Per structure: {relax_times_arr.mean():.3f}s avg, "
        f"{relax_times_arr.std():.3f}s std\n"
        f"  Device: {device}\n"
        f"  Saved to: {out_path}\n{bar}"
    )


if __name__ == "__main__":
    main()
