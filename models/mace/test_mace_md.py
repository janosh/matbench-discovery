"""Run NVT molecular dynamics with a MACE foundation model.

Example:
    uv run python models/mace/test_mace_md.py \
      --input-dir /path/to/reference-trajectories \
      --n-steps 80000
"""

# /// script
# requires-python = ">=3.11,<3.15"
# dependencies = [
#   "ase>=3.25",
#   "mace-torch>=0.3.6",
#   "matbench-discovery>=1.3.1",
#   "pandas>=2.2.2",
#   "torch>=2.5",
#   "tqdm>=4.67",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///

from __future__ import annotations

import argparse
import json
import os
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

from matbench_discovery import today
from matbench_discovery.md import (
    IpiRunConfig,
    MdRunConfig,
    run_md_for_model,
    write_ipi_molecular_crystal_workflows,
)
from matbench_discovery.metrics.md import (
    MdPressureConfig,
    MdRdfConfig,
    MdVdosConfig,
    add_md_combined_error,
    evaluate_pressure_dataset,
    evaluate_rdf_dataset,
    evaluate_reference_dataset,
    evaluate_vdos_dataset,
    write_md_prediction_bundle,
    write_metrics_to_yaml,
)

if TYPE_CHECKING:
    import pandas as pd

SOURCE_IPI_THERMOSTAT_TAU_FS_BY_SYSTEM = {
    "pentacene_295K_Sharma_S": 500.0,
    "picene_295K_Sharma_S": 500.0,
}
SOURCE_IPI_MOLECULAR_CRYSTAL_SYSTEMS = (
    "anthracene_293K_Sharma_S",
    "naphthalene_295K_Sharma_S",
    "pentacene_295K_Sharma_S",
    "picene_295K_Sharma_S",
    "tetracene_295K_Sharma_S",
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing ASE-readable input trajectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for MD trajectories, logs and summary CSV. Defaults to "
            "md-results/<model-key>/<today>-md-nvt."
        ),
    )
    parser.add_argument(
        "--model-key",
        default=os.getenv("MODEL_KEY", "mace-mp-0"),
        help="Model key used in output file names and frame metadata.",
    )
    parser.add_argument(
        "--mace-model",
        default=os.getenv("MACE_MODEL", "medium"),
        help="MACE model name, checkpoint path or checkpoint URL passed to mace_mp().",
    )
    parser.add_argument(
        "--device",
        default=None,
        help=(
            "Device passed to the MACE calculator. Defaults to cuda for generated "
            "i-PI inputs, or auto-detects cuda/cpu for direct MD runs."
        ),
    )
    parser.add_argument(
        "--dtype",
        default="float64",
        choices=("float32", "float64"),
        help="Default floating-point precision for MACE.",
    )
    parser.add_argument("--n-steps", type=int, default=80_000)
    parser.add_argument("--timestep-fs", type=float, default=0.25)
    parser.add_argument("--record-interval", type=int, default=10)
    parser.add_argument("--log-interval", type=int, default=100)
    parser.add_argument(
        "--tdamp-fs",
        type=float,
        default=None,
        help="Nose-Hoover thermostat damping time in fs. Defaults to 100 * timestep.",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing nvt_<model-key>.extxyz files.",
    )
    parser.add_argument(
        "--skip-reference-eval",
        action="store_true",
        help="Skip energy/force RMSE evaluation on the reference MD frames.",
    )
    parser.add_argument(
        "--skip-md",
        action="store_true",
        help="Skip the NVT rollout and only evaluate reference-frame metrics.",
    )
    parser.add_argument(
        "--write-ipi-inputs",
        action="store_true",
        help=(
            "Write molecular-crystal i-PI run folders with run-ase.py, init.xyz, "
            "input.xml and submit.sh, then exit."
        ),
    )
    parser.add_argument(
        "--ipi-output-dir",
        default=None,
        help=(
            "Directory for generated i-PI folders. Defaults to "
            "md-results/<model-key>/<today>-md-ipi."
        ),
    )
    parser.add_argument(
        "--ipi-system",
        action="append",
        default=None,
        help=(
            "System directory name to generate. Can be repeated. Defaults to the "
            "five molecular crystals from the source i-PI workflow."
        ),
    )
    parser.add_argument("--ipi-total-steps", type=int, default=40_000)
    parser.add_argument("--ipi-total-time", type=float, default=40_000)
    parser.add_argument("--ipi-timestep-fs", type=float, default=0.5)
    parser.add_argument("--ipi-output-stride", type=int, default=1)
    parser.add_argument("--ipi-checkpoint-stride", type=int, default=100)
    parser.add_argument("--ipi-thermostat-tau-fs", type=float, default=10)
    parser.add_argument(
        "--no-ipi-source-tau-overrides",
        action="store_true",
        help=(
            "Disable the source i-PI tau overrides for pentacene and picene. "
            "Without this flag, those systems use 500 fs as in the source workflow."
        ),
    )
    parser.add_argument("--ipi-seed", type=int, default=31415)
    parser.add_argument("--ipi-socket-name", default="driver")
    parser.add_argument(
        "--ipi-activate-command",
        default=os.getenv("IPI_ACTIVATE_COMMAND"),
        help="Optional shell command written near the top of each generated submit.sh.",
    )
    parser.add_argument(
        "--ipi-command-prefix",
        default=os.getenv("IPI_COMMAND_PREFIX", "uv run"),
        help="Command prefix for i-pi and Python calls in generated submit.sh.",
    )
    parser.add_argument(
        "--ipi-conda-env",
        default=os.getenv("IPI_CONDA_ENV"),
        help=(
            "Fallback conda environment if --ipi-activate-command is unset. "
            "Kept for compatibility with older generated scripts."
        ),
    )
    parser.add_argument(
        "--cuda-visible-devices",
        default="1",
        help="CUDA_VISIBLE_DEVICES value written to each generated submit.sh.",
    )
    parser.add_argument(
        "--no-unique-ipi-sockets",
        action="store_true",
        help="Use the same i-PI socket name for every generated system.",
    )
    parser.add_argument("--ipi-slurm-job-suffix", default="ipi")
    parser.add_argument("--ipi-nodes", type=int, default=1)
    parser.add_argument("--ipi-cpus-per-task", type=int, default=1)
    parser.add_argument("--ipi-ntasks-per-node", type=int, default=4)
    parser.add_argument("--ipi-walltime", default="48:00:00")
    parser.add_argument("--ipi-partition", default=os.getenv("IPI_PARTITION", "GPU"))
    parser.add_argument("--ipi-gres", default=os.getenv("IPI_GRES", "gpu:v100:1"))
    parser.add_argument(
        "--ipi-tmpdir", default=os.getenv("IPI_TMPDIR", "/home/mjgawkowski/tmp/")
    )
    parser.add_argument("--ipi-sleep-seconds", type=int, default=60)
    parser.add_argument(
        "--write-yaml",
        action="store_true",
        help="Write reference-frame MD metrics to the model YAML file.",
    )
    parser.add_argument(
        "--pred-file-url",
        default=None,
        help="Public URL for the MD reference prediction CSV when writing YAML.",
    )
    parser.add_argument(
        "--run-rdf",
        action="store_true",
        help="Evaluate RDF similarity using MLIP-predicted trajectories.",
    )
    parser.add_argument(
        "--mlip-dir",
        default=None,
        help="Directory containing MLIP nvt_<model-key>.extxyz trajectories. "
        "Defaults to --output-dir.",
    )
    parser.add_argument(
        "--rdf-mode",
        choices=("same", "different"),
        default="same",
        help=(
            "Use matched same-duration trajectories or full different-length "
            "trajectories."
        ),
    )
    parser.add_argument("--rdf-nbins", type=int, default=500)
    parser.add_argument("--ref-frame-dt-fs", type=float, default=None)
    parser.add_argument(
        "--mlip-frame-dt-fs",
        type=float,
        default=None,
        help="MLIP frame spacing in fs. Defaults to timestep * record interval.",
    )
    parser.add_argument("--ref-rdf-settings", default=None)
    parser.add_argument("--mlip-rdf-settings", default=None)
    parser.add_argument(
        "--no-write-rdf-curves",
        action="store_true",
        help="Do not write per-system RDF curve CSVs.",
    )
    parser.add_argument(
        "--run-vdos",
        action="store_true",
        help="Evaluate source-compatible VDOS similarity using MLIP trajectories.",
    )
    parser.add_argument(
        "--vdos-mode",
        choices=("same", "different"),
        default="same",
        help="Use equal post-stride frame counts or full trajectories for VDOS.",
    )
    parser.add_argument(
        "--ref-vdos-frame-dt-fs",
        type=float,
        default=None,
        help="Native reference saved-frame spacing before VDOS stride, in fs.",
    )
    parser.add_argument(
        "--mlip-vdos-frame-dt-fs",
        type=float,
        default=None,
        help="Native MLIP saved-frame spacing before VDOS stride, in fs.",
    )
    parser.add_argument("--ref-vdos-stride", type=int, default=None)
    parser.add_argument("--mlip-vdos-stride", type=int, default=None)
    parser.add_argument("--ref-vdos-padding", type=int, default=None)
    parser.add_argument("--mlip-vdos-padding", type=int, default=None)
    parser.add_argument("--ref-vdos-settings", default=None)
    parser.add_argument("--mlip-vdos-settings", default=None)
    parser.add_argument("--vdos-e-min", type=float, default=None)
    parser.add_argument("--vdos-e-max", type=float, default=None)
    parser.add_argument(
        "--no-normalize-vdos",
        action="store_true",
        help="Disable area normalization before computing VDOS similarity.",
    )
    parser.add_argument(
        "--no-write-vdos-curves",
        action="store_true",
        help="Do not write per-system VDOS spectrum files.",
    )
    parser.add_argument(
        "--run-pressure",
        action="store_true",
        help="Evaluate matched-length pressure histogram percentage error.",
    )
    parser.add_argument("--pressure-bins", type=int, default=80)
    parser.add_argument(
        "--ref-pressure-frame-dt-fs",
        type=float,
        default=None,
        help="Reference saved-frame spacing for matched pressure evaluation, in fs.",
    )
    parser.add_argument(
        "--mlip-pressure-frame-dt-fs",
        type=float,
        default=None,
        help="MLIP saved-frame spacing for matched pressure evaluation, in fs.",
    )
    parser.add_argument("--ref-pressure-settings", default=None)
    parser.add_argument("--mlip-pressure-settings", default=None)
    parser.add_argument(
        "--no-write-pressure-values",
        action="store_true",
        help="Do not write per-system per-frame pressure CSVs.",
    )
    return parser.parse_args()


def main() -> None:
    """Run MACE MD for every trajectory in the input directory."""
    args = parse_args()
    config = MdRunConfig(
        n_steps=args.n_steps,
        timestep_fs=args.timestep_fs,
        record_interval=args.record_interval,
        log_interval=args.log_interval,
        tdamp_fs=args.tdamp_fs,
        seed=args.seed,
        overwrite=args.overwrite,
    )

    if args.write_ipi_inputs:
        ipi_device = args.device or os.getenv("MACE_DEVICE", "cuda")
        ipi_config = IpiRunConfig(
            total_steps=args.ipi_total_steps,
            total_time=args.ipi_total_time,
            timestep_fs=args.ipi_timestep_fs,
            output_stride=args.ipi_output_stride,
            checkpoint_stride=args.ipi_checkpoint_stride,
            thermostat_tau_fs=args.ipi_thermostat_tau_fs,
            thermostat_tau_fs_by_system=(
                None
                if args.no_ipi_source_tau_overrides
                else SOURCE_IPI_THERMOSTAT_TAU_FS_BY_SYSTEM
            ),
            seed=args.ipi_seed,
            socket_name=args.ipi_socket_name,
            mace_model=args.mace_model,
            device=ipi_device,
            dtype=args.dtype,
            dispersion=False,
            command_prefix=args.ipi_command_prefix,
            activate_command=args.ipi_activate_command,
            conda_env=args.ipi_conda_env,
            cuda_visible_devices=args.cuda_visible_devices,
            unique_socket_names=not args.no_unique_ipi_sockets,
            slurm_job_suffix=args.ipi_slurm_job_suffix,
            nodes=args.ipi_nodes,
            cpus_per_task=args.ipi_cpus_per_task,
            ntasks_per_node=args.ipi_ntasks_per_node,
            walltime=args.ipi_walltime,
            partition=args.ipi_partition,
            gres=args.ipi_gres,
            tmpdir=args.ipi_tmpdir,
            sleep_seconds=args.ipi_sleep_seconds,
            overwrite=args.overwrite,
        )
        ipi_output_dir = (
            Path(args.ipi_output_dir)
            if args.ipi_output_dir
            else Path("md-results") / args.model_key / f"{today}-md-ipi"
        )
        df_ipi = write_ipi_molecular_crystal_workflows(
            model_key=args.model_key,
            input_dir=args.input_dir,
            output_dir=ipi_output_dir,
            config=ipi_config,
            md_config=config,
            system_names=args.ipi_system or SOURCE_IPI_MOLECULAR_CRYSTAL_SYSTEMS,
        )
        if df_ipi.empty:
            print("No matching input trajectories found for i-PI generation")
        else:
            print(df_ipi.value_counts("status", dropna=False).to_string())
        print(f"i-PI inputs written to {ipi_output_dir}")
        return

    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else Path("md-results") / args.model_key / f"{today}-md-nvt"
    )
    args.output_dir = str(output_dir)

    import torch
    from mace.calculators import mace_mp

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")

    print(f"Loading MACE model {args.mace_model!r} on {device}")
    calculator = mace_mp(
        model=args.mace_model,
        device=device,
        default_dtype=args.dtype,
    )

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    run_metadata = {
        "model_key": args.model_key,
        "mace_model": args.mace_model,
        "device": device,
        "dtype": args.dtype,
        "input_dir": args.input_dir,
        "output_dir": args.output_dir,
        "versions": {
            dep: version(dep)
            for dep in ("ase", "mace-torch", "matbench-discovery", "torch")
        },
    }
    with open(
        Path(args.output_dir) / f"mace_md_run_{args.model_key}.json",
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(run_metadata, file, indent=2)

    all_md_metrics: dict[str, float] = {}
    md_metric_frames: dict[str, pd.DataFrame] = {}

    if not args.skip_reference_eval:
        _df_ref, md_metrics = evaluate_reference_dataset(
            model_key=args.model_key,
            calculator=calculator,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            config=config,
        )
        md_metric_frames["reference"] = _df_ref
        all_md_metrics |= md_metrics
        print("Reference MD metrics:")
        for key, val in md_metrics.items():
            print(f"  {key}: {val}")

    if not args.skip_md:
        df_results = run_md_for_model(
            model_key=args.model_key,
            calculator=calculator,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            config=config,
        )
        print(df_results.value_counts("status", dropna=False).to_string())

    if args.run_rdf:
        rdf_config = MdRdfConfig(
            mode=args.rdf_mode,
            nbins=args.rdf_nbins,
            ref_frame_dt_fs=args.ref_frame_dt_fs,
            mlip_frame_dt_fs=(
                args.mlip_frame_dt_fs
                if args.mlip_frame_dt_fs is not None
                else args.timestep_fs * args.record_interval
            ),
            ref_settings_path=args.ref_rdf_settings,
            mlip_settings_path=args.mlip_rdf_settings,
            write_rdf_curves=not args.no_write_rdf_curves,
        )
        _df_rdf, rdf_metrics = evaluate_rdf_dataset(
            model_key=args.model_key,
            ref_dir=args.input_dir,
            mlip_dir=args.mlip_dir or args.output_dir,
            output_dir=args.output_dir,
            config=rdf_config,
            md_config=config,
        )
        md_metric_frames["rdf"] = _df_rdf
        all_md_metrics |= rdf_metrics
        print("RDF MD metrics:")
        for key, val in rdf_metrics.items():
            print(f"  {key}: {val}")

    if args.run_vdos:
        vdos_config = MdVdosConfig(
            mode=args.vdos_mode,
            ref_frame_dt_fs=args.ref_vdos_frame_dt_fs,
            mlip_frame_dt_fs=(
                args.mlip_vdos_frame_dt_fs
                if args.mlip_vdos_frame_dt_fs is not None
                else (
                    None
                    if args.mlip_vdos_settings is not None
                    else args.timestep_fs * args.record_interval
                )
            ),
            ref_stride=args.ref_vdos_stride,
            mlip_stride=args.mlip_vdos_stride,
            ref_padding=args.ref_vdos_padding,
            mlip_padding=args.mlip_vdos_padding,
            ref_settings_path=args.ref_vdos_settings,
            mlip_settings_path=args.mlip_vdos_settings,
            e_min=args.vdos_e_min,
            e_max=args.vdos_e_max,
            normalize_area=not args.no_normalize_vdos,
            write_vdos_curves=not args.no_write_vdos_curves,
        )
        _df_vdos, vdos_metrics = evaluate_vdos_dataset(
            model_key=args.model_key,
            ref_dir=args.input_dir,
            mlip_dir=args.mlip_dir or args.output_dir,
            output_dir=args.output_dir,
            config=vdos_config,
            md_config=config,
        )
        md_metric_frames["vdos"] = _df_vdos
        all_md_metrics |= vdos_metrics
        print("VDOS MD metrics:")
        for key, val in vdos_metrics.items():
            print(f"  {key}: {val}")

    if args.run_pressure:
        pressure_config = MdPressureConfig(
            bins=args.pressure_bins,
            ref_frame_dt_fs=args.ref_pressure_frame_dt_fs,
            mlip_frame_dt_fs=(
                args.mlip_pressure_frame_dt_fs
                if args.mlip_pressure_frame_dt_fs is not None
                else (
                    None
                    if args.mlip_pressure_settings is not None
                    else args.timestep_fs * args.record_interval
                )
            ),
            ref_settings_path=args.ref_pressure_settings,
            mlip_settings_path=args.mlip_pressure_settings,
            write_pressure_values=not args.no_write_pressure_values,
        )
        _df_pressure, pressure_metrics = evaluate_pressure_dataset(
            model_key=args.model_key,
            calculator=calculator,
            ref_dir=args.input_dir,
            mlip_dir=args.mlip_dir or args.output_dir,
            output_dir=args.output_dir,
            config=pressure_config,
            md_config=config,
        )
        md_metric_frames["pressure"] = _df_pressure
        all_md_metrics |= pressure_metrics
        print("Pressure MD metrics:")
        for key, val in pressure_metrics.items():
            print(f"  {key}: {val}")

    all_md_metrics = add_md_combined_error(all_md_metrics)
    if "combined_error" in all_md_metrics:
        print(f"Combined MD error: {all_md_metrics['combined_error']}")

    bundled_pred_file_path = None
    if len(md_metric_frames) > 1:
        bundled_pred_file_path = write_md_prediction_bundle(
            model_key=args.model_key,
            output_dir=args.output_dir,
            metric_frames=md_metric_frames,
        )
        print(f"Bundled MD prediction artifact: {bundled_pred_file_path}")

    if args.write_yaml and all_md_metrics:
        from matbench_discovery.enums import Model

        if bundled_pred_file_path is not None:
            pred_file_path = bundled_pred_file_path
        elif args.skip_reference_eval and args.run_rdf:
            mode_name = (
                "same_simulation_length" if args.rdf_mode == "same" else "different"
            )
            pred_file_path = (
                Path(args.output_dir) / f"md_rdf_{mode_name}_{args.model_key}.csv"
            )
        elif args.skip_reference_eval and args.run_vdos:
            mode_name = (
                "same_simulation_length"
                if args.vdos_mode == "same"
                else "different_simulation_length"
            )
            pred_file_path = (
                Path(args.output_dir) / f"md_vdos_{mode_name}_{args.model_key}.csv"
            )
        elif args.skip_reference_eval and args.run_pressure:
            pred_file_path = (
                Path(args.output_dir)
                / f"md_pressure_same_simulation_length_{args.model_key}.csv"
            )
        else:
            pred_file_path = (
                Path(args.output_dir) / f"md_reference_{args.model_key}.csv"
            )
        write_metrics_to_yaml(
            Model(args.model_key),
            all_md_metrics,
            pred_file_path=pred_file_path if args.pred_file_url else None,
            pred_file_url=args.pred_file_url,
        )
        if not args.pred_file_url:
            print("Omitted pred_file from YAML because --pred-file-url is unset")
        print(f"Updated YAML for {args.model_key}")
    print(f"Results written to {args.output_dir}")


if __name__ == "__main__":
    main()
