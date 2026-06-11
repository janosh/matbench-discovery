"""Shared molecular dynamics helpers for model evaluation scripts."""

from __future__ import annotations

import json
import re
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING, Any

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.md import MDLogger
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from tqdm import tqdm

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase.calculators.calculator import Calculator


@dataclass(frozen=True)
class MdRunConfig:
    """Configuration for NVT molecular dynamics runs."""

    n_steps: int = 80000
    timestep_fs: float = 0.25
    record_interval: int = 10
    log_interval: int = 100
    tdamp_fs: float | None = None
    seed: int | None = None
    overwrite: bool = False
    skip_dir_prefixes: tuple[str, ...] = ("H_",)
    skip_filenames: tuple[str, ...] = (
        "isolated_atom_C.extxyz",
        "isolated_atom_C.xyz",
        "isolated_atom_H.extxyz",
        "isolated_atom_H.xyz",
    )


@dataclass(frozen=True)
class MdRunResult:
    """File paths and metadata produced by one MD run."""

    model_key: str
    input_file: str
    system_name: str
    temperature_k: float
    traj_file: str
    log_file: str
    n_steps: int
    timestep_fs: float
    record_interval: int
    frames_written: int
    status: str = "completed"
    error: str | None = None


@dataclass(frozen=True)
class IpiRunConfig:
    """Configuration for generating i-PI molecular-crystal run folders."""

    total_steps: int = 40_000
    total_time: float = 40_000
    timestep_fs: float = 0.5
    output_stride: int = 1
    checkpoint_stride: int = 100
    thermostat_tau_fs: float = 10
    thermostat_tau_fs_by_system: dict[str, float] | None = None
    seed: int = 31415
    socket_name: str = "driver"
    mace_model: str = "medium"
    device: str = "cuda"
    dtype: str = "float64"
    dispersion: bool = False
    command_prefix: str | None = "uv run"
    activate_command: str | None = None
    conda_env: str | None = None
    cuda_visible_devices: str | None = "1"
    unique_socket_names: bool = True
    slurm_job_suffix: str = "ipi"
    nodes: int = 1
    cpus_per_task: int = 1
    ntasks_per_node: int = 4
    walltime: str = "48:00:00"
    partition: str = "GPU"
    gres: str = "gpu:v100:1"
    tmpdir: str | None = "/home/mjgawkowski/tmp/"
    sleep_seconds: int = 60
    overwrite: bool = False


@dataclass(frozen=True)
class IpiRunFiles:
    """File paths and metadata produced for one i-PI run folder."""

    model_key: str
    input_file: str
    system_name: str
    temperature_k: float
    run_dir: str
    init_file: str
    input_xml_file: str
    run_ase_file: str
    submit_file: str
    status: str = "completed"
    error: str | None = None


def find_trajectory_files(
    input_dir: str | Path,
    *,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> list[Path]:
    """Find trajectory files in a directory tree.

    Args:
        input_dir: Directory to scan recursively.
        extensions: File extensions to include.

    Returns:
        Sorted list of matching file paths.
    """
    root = Path(input_dir)
    if not root.is_dir():
        raise NotADirectoryError(f"{root} is not a directory")

    normalized_extensions = tuple(
        ext if ext.startswith(".") else f".{ext}" for ext in extensions
    )
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file()
        and any(path.name.endswith(ext) for ext in normalized_extensions)
    )


def read_trajectory(file_path: str | Path) -> list[Atoms]:
    """Read all frames from an ASE-readable trajectory file."""
    frames = ase.io.read(file_path, index=":")
    if isinstance(frames, Atoms):
        return [frames]
    return list(frames)


def extract_temperature_from_name(name: str) -> float | None:
    """Extract a temperature in Kelvin from strings like ``CsPbBr3_300K``."""
    if match := re.search(r"(?P<temperature>\d+(?:\.\d+)?)K(?![A-Za-z0-9.])", name):
        return float(match["temperature"])
    return None


def infer_temperature(
    file_path: str | Path, atoms: Atoms | None = None
) -> float | None:
    """Infer MD temperature from Atoms metadata or path components."""
    if atoms is not None:
        for key in ("temperature_K", "temperature", "tempK", "temp"):
            if key in atoms.info:
                try:
                    return float(atoms.info[key])
                except (TypeError, ValueError):
                    continue

    path = Path(file_path)
    for part in (path.stem, *[parent.name for parent in path.parents]):
        temperature = extract_temperature_from_name(part)
        if temperature is not None:
            return temperature
    return None


def should_skip_trajectory(file_path: str | Path, config: MdRunConfig) -> bool:
    """Return whether a trajectory should be skipped for this MD task."""
    path = Path(file_path)
    return path.name in config.skip_filenames or any(
        path.parent.name.startswith(prefix) for prefix in config.skip_dir_prefixes
    )


def _format_ipi_number(value: float) -> str:
    """Format XML scalar values without noisy trailing zeros."""
    return f"{value:g}"


def _ipi_input_preference(file_path: Path) -> tuple[int, str]:
    """Rank candidate initial structures when one system has multiple files."""
    preferred_names = {
        "traj.extxyz": 0,
        "traj.xyz": 1,
        "out.extxyz": 2,
        "out.xyz": 3,
    }
    if file_path.name in preferred_names:
        return preferred_names[file_path.name], file_path.name
    if file_path.name.endswith("_trajectory.extxyz"):
        return 4, file_path.name
    if file_path.name.endswith("_trajectory.xyz"):
        return 5, file_path.name
    if file_path.name.endswith(".extxyz"):
        return 6, file_path.name
    if file_path.name.endswith(".xyz"):
        return 7, file_path.name
    return 8, file_path.name


def _safe_shell_name(name: str) -> str:
    """Return a shell/socket-safe identifier fragment."""
    return re.sub(r"[^A-Za-z0-9-]+", "-", name).strip("-").lower()


def _one_ipi_input_per_system(file_paths: Sequence[Path]) -> list[Path]:
    """Keep the preferred initial-structure file for each system directory."""
    best_by_system: dict[str, Path] = {}
    for file_path in file_paths:
        system_name = file_path.parent.name
        current = best_by_system.get(system_name)
        if current is None or _ipi_input_preference(file_path) < _ipi_input_preference(
            current
        ):
            best_by_system[system_name] = file_path
    return sorted(best_by_system.values())


def _ipi_socket_name(system_name: str, config: IpiRunConfig) -> str:
    """Return the Unix socket name for one generated i-PI run."""
    if not config.unique_socket_names:
        return config.socket_name
    molecule = system_name.split("_", maxsplit=1)[0]
    safe_molecule = _safe_shell_name(molecule).replace("-", "_")
    safe_socket_base = _safe_shell_name(config.socket_name).replace("-", "_")
    return f"{safe_socket_base}_{safe_molecule}"[:96]


def _thermostat_tau_fs(system_name: str, config: IpiRunConfig) -> float:
    """Return the configured thermostat tau for one molecular crystal."""
    if config.thermostat_tau_fs_by_system:
        molecule = system_name.split("_", maxsplit=1)[0]
        for key in (system_name, molecule):
            if key in config.thermostat_tau_fs_by_system:
                return config.thermostat_tau_fs_by_system[key]
    return config.thermostat_tau_fs


def _render_ipi_input_xml(
    system_name: str, temperature_k: float, config: IpiRunConfig, socket_name: str
) -> str:
    """Render an i-PI input file matching the molecular-crystal workflow."""
    temperature = _format_ipi_number(temperature_k)
    thermostat_tau = _format_ipi_number(_thermostat_tau_fs(system_name, config))
    properties = (
        "[ step, time{picosecond}, conserved{electronvolt}, "
        "temperature{kelvin}, kinetic_md{electronvolt}, "
        "potential{electronvolt}, pressure_md{megapascal} ]"
    )
    return dedent(
        f"""\
        <simulation mode='static' verbosity='medium'>
          <output prefix='simulation'>
            <properties filename='out' stride='{config.output_stride}'>
              {properties}
            </properties>
            <trajectory filename='pos' stride='{config.output_stride}' format='ase'>
              positions
            </trajectory>
            <trajectory filename='vel' stride='{config.output_stride}' format='xyz'>
              velocities
            </trajectory>
            <trajectory filename='for' stride='{config.output_stride}' format='xyz'>
              forces
            </trajectory>
            <checkpoint stride='{config.checkpoint_stride}'/>
          </output>
          <total_time> {_format_ipi_number(config.total_time)} </total_time>
          <total_steps>{config.total_steps}</total_steps>
          <prng>
            <seed>{config.seed}</seed>
          </prng>
          <ffsocket name='{socket_name}' mode='unix'>
            <address> {socket_name} </address>
          </ffsocket>
          <system>
            <initialize nbeads='1'>
              <file mode='ase'> init.xyz </file>
              <velocities mode='thermal' units='kelvin'> {temperature} </velocities>
            </initialize>
            <forces>
              <force forcefield='{socket_name}'> </force>
            </forces>
            <motion mode='dynamics'>
              <dynamics mode='nvt'>
                <thermostat mode='langevin'>
                  <tau units='femtosecond'>
                    {thermostat_tau}
                  </tau>
                </thermostat>
                <timestep units='femtosecond'>
                  {_format_ipi_number(config.timestep_fs)}
                </timestep>
              </dynamics>
            </motion>
            <ensemble>
              <temperature units='kelvin'> {temperature} </temperature>
            </ensemble>
          </system>
        </simulation>
        """
    )


def _render_mace_ipi_runner(config: IpiRunConfig, socket_name: str) -> str:
    """Render the ASE socket client used by i-PI to call MACE-MP."""
    dispersion = "true" if config.dispersion else "false"
    return dedent(
        f'''\
        """ASE socket client for i-PI + MACE-MP."""

        from __future__ import annotations

        import os

        from ase.calculators.socketio import SocketClient
        from ase.io import read
        from mace.calculators import mace_mp


        def env_bool(name: str, default: str) -> bool:
            """Read a bool-like environment variable."""
            return os.getenv(name, default).lower() in {{"1", "true", "yes", "on"}}


        atoms = read("init.xyz", index=0)
        atoms.calc = mace_mp(
            model=os.getenv("MACE_MODEL", "{config.mace_model}"),
            device=os.getenv("MACE_DEVICE", "{config.device}"),
            dispersion=env_bool("MACE_DISPERSION", "{dispersion}"),
            default_dtype=os.getenv("MACE_DTYPE", "{config.dtype}"),
        )

        client = SocketClient(
            unixsocket=os.getenv("IPI_SOCKET", "{socket_name}")
        )
        client.run(atoms, use_stress=True)
        '''
    )


def _render_ipi_submit_script(config: IpiRunConfig, socket_name: str) -> str:
    """Render a plain bash launcher for one i-PI molecular-crystal run."""
    activate_command = (
        config.activate_command
        if config.activate_command is not None
        else (f"conda activate {config.conda_env}" if config.conda_env else "")
    )
    activate_block = f"\n{activate_command}\n" if activate_command else ""
    command_prefix = f"{config.command_prefix} " if config.command_prefix else ""
    tmpdir_export = f"export TMPDIR={config.tmpdir}" if config.tmpdir else ""
    cuda_export = (
        f"export CUDA_VISIBLE_DEVICES={config.cuda_visible_devices}"
        if config.cuda_visible_devices is not None
        else "# CUDA_VISIBLE_DEVICES left unchanged"
    )
    return dedent(
        f"""\
        #!/bin/bash -l
        set -euo pipefail
        {activate_block}
        {tmpdir_export}
        export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-1}}
        export CRAY_CUDA_MPS=1
        export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
        {cuda_export}
        export IPI_SOCKET={socket_name}
        export MACE_MODEL={config.mace_model}
        export MACE_DEVICE={config.device}
        export MACE_DTYPE={config.dtype}
        export MACE_DISPERSION={"true" if config.dispersion else "false"}

        rm -f /tmp/ipi_${{IPI_SOCKET}}*

        if [ -e RESTART ]; then
            echo "RESTART file found. Running i-PI with RESTART."
            {command_prefix}i-pi RESTART &> log.ipi &
        else
            echo "RESTART file not found. Running i-PI with input.xml."
            {command_prefix}i-pi input.xml &> log.ipi &
        fi

        sleep {config.sleep_seconds}

        {command_prefix}python3 run-ase.py &

        wait
        """
    )


def _render_ipi_submit_all_script(model_key: str) -> str:
    """Render a root-level launcher for generated i-PI submit scripts."""
    return dedent(
        f"""\
        #!/bin/bash -l
        set -euo pipefail

        MODE="${{1:-parallel}}"
        ROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"

        mapfile -t scripts < <(find "$ROOT" -path '*/{model_key}/submit.sh' | sort)
        if [ "${{#scripts[@]}}" -eq 0 ]; then
            echo "No submit.sh files found under $ROOT"
            exit 1
        fi

        for script in "${{scripts[@]}}"; do
            run_dir="$(dirname "$script")"
            echo "Launching $run_dir"
            if [ "$MODE" = "sequential" ]; then
                (cd "$run_dir" && bash submit.sh)
            else
                (
                    cd "$run_dir"
                    nohup bash submit.sh > submit.log 2>&1 &
                    echo "$!" > submit.pid
                )
            fi
        done

        if [ "$MODE" != "sequential" ]; then
            echo "Launched ${{#scripts[@]}} run(s) in the background."
            echo "Logs: */{model_key}/submit.log"
            echo "PIDs: */{model_key}/submit.pid"
        fi
        """
    )


def write_ipi_molecular_crystal_workflow(
    *,
    model_key: str,
    input_file: str | Path,
    init_structure: Atoms,
    temperature_k: float,
    output_dir: str | Path,
    config: IpiRunConfig,
) -> IpiRunFiles:
    """Write ``run-ase.py``, ``init.xyz``, ``input.xml`` and ``submit.sh``."""
    system_name = Path(input_file).parent.name
    run_dir = Path(output_dir) / system_name / model_key
    run_dir.mkdir(parents=True, exist_ok=True)

    init_path = run_dir / "init.xyz"
    input_xml_path = run_dir / "input.xml"
    run_ase_path = run_dir / "run-ase.py"
    submit_path = run_dir / "submit.sh"
    target_files = (init_path, input_xml_path, run_ase_path, submit_path)
    socket_name = _ipi_socket_name(system_name, config)

    if any(path.exists() for path in target_files) and not config.overwrite:
        return IpiRunFiles(
            model_key=model_key,
            input_file=str(input_file),
            system_name=system_name,
            temperature_k=temperature_k,
            run_dir=str(run_dir),
            init_file=str(init_path),
            input_xml_file=str(input_xml_path),
            run_ase_file=str(run_ase_path),
            submit_file=str(submit_path),
            status="skipped_existing",
        )

    ase.io.write(init_path, init_structure, format="extxyz")
    input_xml_path.write_text(
        _render_ipi_input_xml(system_name, temperature_k, config, socket_name),
        encoding="utf-8",
    )
    run_ase_path.write_text(
        _render_mace_ipi_runner(config, socket_name), encoding="utf-8"
    )
    submit_path.write_text(
        _render_ipi_submit_script(config, socket_name), encoding="utf-8"
    )
    submit_path.chmod(0o755)

    return IpiRunFiles(
        model_key=model_key,
        input_file=str(input_file),
        system_name=system_name,
        temperature_k=temperature_k,
        run_dir=str(run_dir),
        init_file=str(init_path),
        input_xml_file=str(input_xml_path),
        run_ase_file=str(run_ase_path),
        submit_file=str(submit_path),
    )


def write_ipi_molecular_crystal_workflows(
    *,
    model_key: str,
    input_dir: str | Path,
    output_dir: str | Path,
    config: IpiRunConfig,
    md_config: MdRunConfig | None = None,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
    system_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Generate i-PI run folders for every molecular-crystal input trajectory."""
    md_config = md_config or MdRunConfig()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    run_params = {
        "model_key": model_key,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": asdict(config),
        "extensions": list(extensions),
        "system_names": list(system_names) if system_names is not None else None,
    }
    with open(
        output_path / f"ipi_run_params_{model_key}.json", mode="w", encoding="utf-8"
    ) as file:
        json.dump(run_params, file, indent=2)

    results: list[dict[str, Any]] = []
    trajectory_files = find_trajectory_files(input_dir, extensions=extensions)
    if system_names is not None:
        allowed_system_names = set(system_names)
        trajectory_files = [
            file_path
            for file_path in trajectory_files
            if file_path.parent.name in allowed_system_names
        ]
    trajectory_files = _one_ipi_input_per_system(trajectory_files)
    for file_path in tqdm(
        trajectory_files, desc=f"Writing i-PI inputs for {model_key}"
    ):
        system_name = file_path.parent.name
        if should_skip_trajectory(file_path, md_config):
            results.append(
                {
                    "model_key": model_key,
                    "input_file": str(file_path),
                    "system_name": system_name,
                    "status": "skipped_input",
                }
            )
            continue

        try:
            frames = read_trajectory(file_path)
            if not frames:
                raise ValueError(f"No frames found in {file_path}")

            temperature_k = infer_temperature(file_path, atoms=frames[0])
            if temperature_k is None:
                raise ValueError(f"Could not infer temperature from {file_path}")

            result = write_ipi_molecular_crystal_workflow(
                model_key=model_key,
                input_file=file_path,
                init_structure=frames[0],
                temperature_k=temperature_k,
                output_dir=output_path,
                config=config,
            )
            results.append(asdict(result))
        except Exception as exc:  # noqa: BLE001 - continue writing other inputs
            results.append(
                {
                    "model_key": model_key,
                    "input_file": str(file_path),
                    "system_name": system_name,
                    "status": "failed",
                    "error": repr(exc),
                    "traceback": traceback.format_exc(),
                }
            )

    df_results = pd.DataFrame(results)
    df_results.to_csv(output_path / f"ipi_summary_{model_key}.csv", index=False)
    submit_all_path = output_path / "submit_all.sh"
    submit_all_path.write_text(
        _render_ipi_submit_all_script(model_key), encoding="utf-8"
    )
    submit_all_path.chmod(0o755)
    return df_results


def _init_velocities(atoms: Atoms, temperature_k: float, seed: int | None) -> None:
    kwargs: dict[str, Any] = {}
    if seed is not None:
        kwargs["rng"] = np.random.default_rng(seed)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_k, **kwargs)


def _make_nvt_dynamics(
    atoms: Atoms, temperature_k: float, config: MdRunConfig
) -> NoseHooverChainNVT:
    timestep = config.timestep_fs * units.fs
    tdamp_fs = config.tdamp_fs or 100 * config.timestep_fs
    return NoseHooverChainNVT(
        atoms,
        temperature_K=temperature_k,
        timestep=timestep,
        tdamp=tdamp_fs * units.fs,
        tchain=1,
    )


def run_nvt_simulation(
    *,
    init_structure: Atoms,
    temperature_k: float,
    calculator: Calculator,
    model_key: str,
    system_name: str,
    input_file: str | Path,
    output_dir: str | Path,
    config: MdRunConfig,
) -> MdRunResult:
    """Run one NVT MD simulation and stream frames to an extxyz trajectory."""
    system_dir = Path(output_dir) / system_name
    system_dir.mkdir(parents=True, exist_ok=True)

    traj_path = system_dir / f"nvt_{model_key}.extxyz"
    log_path = system_dir / f"md_{model_key}.log"
    if traj_path.exists():
        if config.overwrite:
            traj_path.unlink()
        else:
            return MdRunResult(
                model_key=model_key,
                input_file=str(input_file),
                system_name=system_name,
                temperature_k=temperature_k,
                traj_file=str(traj_path),
                log_file=str(log_path),
                n_steps=config.n_steps,
                timestep_fs=config.timestep_fs,
                record_interval=config.record_interval,
                frames_written=0,
                status="skipped_existing",
            )

    atoms = init_structure.copy()
    atoms.calc = calculator
    _init_velocities(atoms, temperature_k, config.seed)
    dyn = _make_nvt_dynamics(atoms, temperature_k, config)

    dyn.attach(
        MDLogger(
            dyn,
            atoms,
            str(log_path),
            header=True,
            stress=False,
            peratom=False,
            mode="w",
        ),
        interval=config.log_interval,
    )

    frames_written = 0
    n_frames_expected = config.n_steps // config.record_interval + 1

    with tqdm(
        total=n_frames_expected,
        desc=f"MD frames: {system_name}",
        unit="frame",
        leave=False,
    ) as progress_bar:

        def store_frame() -> None:
            nonlocal frames_written
            frame = atoms.copy()
            frame.calc = None
            frame.info["model_key"] = model_key
            frame.info["temperature_K"] = temperature_k
            frame.info["md_step"] = int(getattr(dyn, "nsteps", frames_written))
            ase.io.write(str(traj_path), frame, append=traj_path.exists())
            frames_written += 1
            progress_bar.update()

        dyn.attach(store_frame, interval=config.record_interval)
        dyn.run(config.n_steps)

    return MdRunResult(
        model_key=model_key,
        input_file=str(input_file),
        system_name=system_name,
        temperature_k=temperature_k,
        traj_file=str(traj_path),
        log_file=str(log_path),
        n_steps=config.n_steps,
        timestep_fs=config.timestep_fs,
        record_interval=config.record_interval,
        frames_written=frames_written,
    )


def run_md_for_model(
    *,
    model_key: str,
    calculator: Calculator,
    input_dir: str | Path,
    output_dir: str | Path,
    config: MdRunConfig,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> pd.DataFrame:
    """Run NVT MD for every trajectory in ``input_dir`` with one model calculator.

    The first frame of each input trajectory is used as the initial structure. A
    ``md_summary_<model_key>.csv`` file and ``run_params_<model_key>.json`` file are
    written to ``output_dir``.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    run_params = {
        "model_key": model_key,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": asdict(config),
        "extensions": list(extensions),
    }
    with open(
        output_path / f"run_params_{model_key}.json", mode="w", encoding="utf-8"
    ) as file:
        json.dump(run_params, file, indent=2)

    trajectory_files = find_trajectory_files(input_dir, extensions=extensions)
    results: list[dict[str, Any]] = []

    for file_path in tqdm(trajectory_files, desc=f"Running MD with {model_key}"):
        system_name = file_path.parent.name
        if should_skip_trajectory(file_path, config):
            results.append(
                {
                    "model_key": model_key,
                    "input_file": str(file_path),
                    "system_name": system_name,
                    "status": "skipped_input",
                }
            )
            continue

        try:
            frames = read_trajectory(file_path)
            if not frames:
                raise ValueError(f"No frames found in {file_path}")

            temperature_k = infer_temperature(file_path, atoms=frames[0])
            if temperature_k is None:
                raise ValueError(f"Could not infer temperature from {file_path}")

            result = run_nvt_simulation(
                init_structure=frames[0],
                temperature_k=temperature_k,
                calculator=calculator,
                model_key=model_key,
                system_name=system_name,
                input_file=file_path,
                output_dir=output_path,
                config=config,
            )
            results.append(asdict(result))
        except Exception as exc:  # noqa: BLE001 - keep batch MD runs fault-tolerant
            results.append(
                {
                    "model_key": model_key,
                    "input_file": str(file_path),
                    "system_name": system_name,
                    "status": "failed",
                    "error": repr(exc),
                    "traceback": traceback.format_exc(),
                }
            )

    df_results = pd.DataFrame(results)
    df_results.to_csv(output_path / f"md_summary_{model_key}.csv", index=False)
    return df_results
