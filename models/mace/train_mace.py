"""This example training script demonstrates how to train a MACE model
using multiple GPUs with PyTorch's DistributedDataParallel.

earlier versions of this script were used to train the
`2023-08-14-mace-yuan-trained-mptrj-04` and
`2023-07-14-mace-universal-2-big-128-6` models.

If you want to fine-tune an existing MACE checkpoint rather than train a
model from scratch, install the foundations branch instead which has an
interface just for that.
pip install git+https://github.com/ACEsuit/mace@foundations

The training here did _not_ use MACE's newest attention block feature
which in our testing performed significantly worse than
`RealAgnosticResidualInteractionBlock`.
"""

import ast
import json
import os
from pathlib import Path
from typing import Any

import mace
import mace.data
import numpy as np
import torch.distributed
import torch.nn.functional
from e3nn import o3
from mace import modules, tools
from mace.mace.data import HDF5Dataset
from mace.tools import torch_geometric
from mace.tools.scripts_utils import (
    LRScheduler,
    create_error_table,
    get_atomic_energies,
    get_config_type_weights,
    get_dataset_from_xyz,
    get_files_with_suffix,
    get_loss_fn,
)
from mace.tools.slurm_distributed import DistributedEnvironment
from torch.nn.parallel import DistributedDataParallel
from torch.optim.swa_utils import SWALR, AveragedModel
from torch_ema import ExponentialMovingAverage

from matbench_discovery import WANDB_PATH, today
from matbench_discovery.slurm import slurm_submit

__author__ = "Yuan Chiang, Ilyes Batatia, Gregor Simm, David Kovacs"
__date__ = "2023-09-18"


module_dir = os.path.dirname(__file__)

slurm_vars = slurm_submit(
    job_name=(job_name := "train-mace-mptrj"),
    account="matgen",
    time="7-00:00:00",
    out_dir=os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}"),
    slurm_flags="""
        -q preempt
        -C gpu
        -G 40
        -N 10
        --ntasks=40
        --ntasks-per-node=4
        --cpus-per-task=8
        --time=02:00:00
        --time-min=01:00:00
        --comment=7-00:00:00
        --signal=B:USR1@180
        --requeue
        --exclusive
        --open-mode=append""",
    pre_cmd=". ~/.venv/py311/bin/activate",
)


def main(**kwargs: Any) -> None:
    """Main function for training MACE models."""
    args = tools.build_default_arg_parser()
    args.parse_args(["--name", kwargs["name"], "--train_file", kwargs["train_file"]])

    for key, values in kwargs.items():
        setattr(args, key, values)

    tag = tools.get_tag(name=args.name, seed=args.seed)
    if args.distributed:
        try:
            distr_env = DistributedEnvironment()
        except Exception as exc:
            print(f"Error specifying environment for distributed training: {exc}")
            return
        world_size = distr_env.world_size
        local_rank = distr_env.local_rank
        rank = distr_env.rank
        if rank == 0:
            print(distr_env)
        torch.distributed.init_process_group(backend="nccl")
    else:
        rank = 0
        local_rank = world_size = 1

    # Setup
    tools.set_seeds(args.seed)
    tools.setup_logger(level=args.log_level, tag=tag, directory=args.log_dir, rank=rank)

    if args.distributed:
        torch.cuda.set_device(local_rank)
        print(f"Process group initialized: {torch.distributed.is_initialized()}")
        print(f"Processes: {world_size}")

    try:
        print(f"MACE version: {mace.__version__}")
    except AttributeError:
        print("Cannot find MACE version, please install MACE via pip")
    print(f"Configuration: {args}")

    tools.set_default_dtype(args.default_dtype)
    device = tools.init_device(args.device)

    if args.statistics_file is not None:
        with open(args.statistics_file) as file:
            statistics = json.load(file)
        print("Using statistics json file")
        args.r_max = statistics["r_max"]
        args.atomic_numbers = statistics["atomic_numbers"]
        args.mean = statistics["mean"]
        args.std = statistics["std"]
        args.avg_num_neighbors = statistics["avg_num_neighbors"]
        args.compute_avg_num_neighbors = False
        args.E0s = statistics["atomic_energies"]

    # Data preparation
    if args.train_file.endswith(".xyz"):
        if args.valid_file is not None and not args.valid_file.endswith(".xyz"):
            raise RuntimeError(
                f"valid_file must be .xyz if train_file is .xyz, got {args.valid_file}"
            )
        config_type_weights = get_config_type_weights(args.config_type_weights)
        collections, atomic_energies_dict = get_dataset_from_xyz(
            train_path=args.train_file,
            valid_path=args.valid_file,
            valid_fraction=args.valid_fraction,
            config_type_weights=config_type_weights,
            test_path=args.test_file,
            seed=args.seed,
            energy_key=args.energy_key,
            forces_key=args.forces_key,
            stress_key=args.stress_key,
            virials_key=args.virials_key,
            dipole_key=args.dipole_key,
            charges_key=args.charges_key,
        )

        test_config_lens = ", ".join(
            f"{name}: {len(test_configs)}" for name, test_configs in collections.tests
        )
        print(
            f"Total number of configurations: train={len(collections.train)}, valid="
            f"{len(collections.valid)}, tests=[{test_config_lens}]"
        )
    elif args.train_file.endswith(".h5"):
        atomic_energies_dict = collections = None
    else:
        raise RuntimeError(
            f"train_file must be either .xyz or .h5, got {args.train_file}"
        )

    # Atomic number table
    if args.atomic_numbers is None:
        if not args.train_file.endswith(".xyz"):
            raise RuntimeError(
                "atomic_numbers must be provided if train_file is not in .xyz format"
            )
        z_table = tools.get_atomic_number_table_from_zs(
            z
            for configs in (collections.train, collections.valid)
            for config in configs
            for z in config.atomic_numbers
        )
    else:
        if args.statistics_file is None:
            print("Using atomic numbers from command line argument")
        else:
            print("Using atomic numbers from statistics file")
        zs_list = ast.literal_eval(args.atomic_numbers)
        if not isinstance(zs_list, list):
            raise ValueError("atomic_numbers did not parse to a list")
        z_table = tools.get_atomic_number_table_from_zs(zs_list)
    print(z_table)

    if atomic_energies_dict is None or len(atomic_energies_dict) == 0:
        if args.train_file.endswith(".xyz"):
            atomic_energies_dict = get_atomic_energies(
                args.E0s, collections.train, z_table
            )
        else:
            atomic_energies_dict = get_atomic_energies(args.E0s, None, z_table)

    if args.model == "AtomicDipolesMACE":
        atomic_energies: np.ndarray | None = None
        dipole_only = True
        compute_dipole = True
        compute_energy = False
        args.compute_forces = False
        compute_virials = False
        args.compute_stress = False
    else:
        dipole_only = False
        if args.model == "EnergyDipolesMACE":
            compute_dipole = True
            compute_energy = True
            args.compute_forces = True
            compute_virials = False
            args.compute_stress = False
        else:
            compute_energy = True
            compute_dipole = False

        atomic_energies = np.array([atomic_energies_dict[z] for z in z_table.zs])
        print(f"Atomic energies: {atomic_energies.tolist()}")

    if args.train_file.endswith(".xyz"):
        train_set = [
            mace.data.AtomicData.from_config(config, z_table=z_table, cutoff=args.r_max)
            for config in collections.train
        ]
        valid_set = [
            mace.data.AtomicData.from_config(config, z_table=z_table, cutoff=args.r_max)
            for config in collections.valid
        ]
    else:
        train_set = HDF5Dataset(args.train_file, r_max=args.r_max, z_table=z_table)
        valid_set = HDF5Dataset(args.valid_file, r_max=args.r_max, z_table=z_table)

    train_sampler, valid_sampler = None, None
    if args.distributed:
        train_sampler = torch.utils.mace.data.distributed.DistributedSampler(
            train_set,
            num_replicas=world_size,
            rank=rank,
            shuffle=True,
            drop_last=True,
            seed=args.seed,
        )
        valid_sampler = torch.utils.data.distributed.DistributedSampler(
            valid_set,
            num_replicas=world_size,
            rank=rank,
            shuffle=True,
            drop_last=True,
            seed=args.seed,
        )

    train_loader = torch_geometric.dataloader.DataLoader(
        dataset=train_set,
        batch_size=args.batch_size,
        sampler=train_sampler,
        shuffle=(train_sampler is None),
        drop_last=False,
        pin_memory=args.pin_memory,
        num_workers=args.num_workers,
    )
    valid_loader = torch_geometric.dataloader.DataLoader(
        dataset=valid_set,
        batch_size=args.valid_batch_size,
        sampler=valid_sampler,
        shuffle=(valid_sampler is None),
        drop_last=False,
        pin_memory=args.pin_memory,
        num_workers=args.num_workers,
    )

    loss_fn: torch.nn.Module = get_loss_fn(
        args.loss,
        args.energy_weight,
        args.forces_weight,
        args.stress_weight,
        args.virials_weight,
        args.dipole_weight,
        dipole_only,
        compute_dipole,
    )
    print(loss_fn)

    if args.compute_avg_num_neighbors:
        args.avg_num_neighbors = modules.compute_avg_num_neighbors(train_loader)
    print(f"Average number of neighbors: {args.avg_num_neighbors}")

    # Selecting outputs
    compute_virials = False
    if args.loss in ("stress", "virials", "huber"):
        compute_virials = True
        args.compute_stress = True
        args.error_table = "PerAtomRMSEstressvirials"

    if args.loss == "uip" and args.compute_stress:
        compute_virials = True

    output_args = {
        "energy": compute_energy,
        "forces": args.compute_forces,
        "virials": compute_virials,
        "stress": args.compute_stress,
        "dipoles": compute_dipole,
    }
    print(f"Selected the following outputs: {output_args}")

    # Build model
    print("Building model")
    if args.num_channels is not None and args.max_L is not None:
        if args.num_channels <= 0:
            raise ValueError("num_channels must be a positive integer")
        if args.max_L < 0:
            raise ValueError("max_L must be a non-negative integer")
        args.hidden_irreps = o3.Irreps(
            (args.num_channels * o3.Irreps.spherical_harmonics(args.max_L))
            .sort()
            .irreps.simplify()
        )

    if len({irrep.mul for irrep in o3.Irreps(args.hidden_irreps)}) != 1:
        raise ValueError(
            "All channels must have the same dimension, use the num_channels and max_L"
            " keywords to specify the number of channels and the maximum L"
        )

    print(f"Hidden irreps: {args.hidden_irreps}")

    model_config = dict(
        r_max=args.r_max,
        num_bessel=args.num_radial_basis,
        num_polynomial_cutoff=args.num_cutoff_basis,
        max_ell=args.max_ell,
        interaction_cls=modules.interaction_classes[args.interaction],
        num_interactions=args.num_interactions,
        num_elements=len(z_table),
        hidden_irreps=o3.Irreps(args.hidden_irreps),
        atomic_energies=atomic_energies,
        avg_num_neighbors=args.avg_num_neighbors,
        atomic_numbers=z_table.zs,
    )

    model: torch.nn.Module

    if args.scaling == "no_scaling":
        args.std = 1.0
        print("No scaling selected")
    elif args.mean is None or args.std is None:
        args.mean, args.std = modules.scaling_classes[args.scaling](
            train_loader, atomic_energies
        )

    if args.model == "MACE":
        model = modules.ScaleShiftMACE(
            **model_config,
            correlation=args.correlation,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[
                "RealAgnosticInteractionBlock"
            ],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
            atomic_inter_scale=args.std,
            atomic_inter_shift=0.0,
            radial_MLP=ast.literal_eval(args.radial_MLP),
            radial_type=args.radial_type,
        )
    elif args.model == "ScaleShiftMACE":
        model = modules.ScaleShiftMACE(
            **model_config,
            correlation=args.correlation,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[args.interaction_first],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
            atomic_inter_scale=args.std,
            atomic_inter_shift=args.mean,
            radial_MLP=ast.literal_eval(args.radial_MLP),
            radial_type=args.radial_type,
        )
    elif args.model == "ScaleShiftBOTNet":
        model = modules.ScaleShiftBOTNet(
            **model_config,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[args.interaction_first],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
            atomic_inter_scale=args.std,
            atomic_inter_shift=args.mean,
        )
    elif args.model == "BOTNet":
        model = modules.BOTNet(
            **model_config,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[args.interaction_first],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
        )
    elif args.model == "AtomicDipolesMACE":
        # std_df = modules.scaling_classes["rms_dipoles_scaling"](train_loader)
        if args.loss != "dipole":
            raise ValueError("Must use dipole loss with AtomicDipolesMACE model")
        if args.error_table != "DipoleRMSE":
            raise ValueError(
                "Must use error_table DipoleRMSE with AtomicDipolesMACE model"
            )
        model = modules.AtomicDipolesMACE(
            **model_config,
            correlation=args.correlation,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[
                "RealAgnosticInteractionBlock"
            ],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
            # dipole_scale=1,
            # dipole_shift=0,
        )
    elif args.model == "EnergyDipolesMACE":
        # std_df = modules.scaling_classes["rms_dipoles_scaling"](train_loader)
        if args.loss != "energy_forces_dipole":
            raise ValueError(
                "Use energy_forces_dipole loss with EnergyDipolesMACE model"
            )
        if args.error_table != "EnergyDipoleRMSE":
            raise ValueError(
                "Use error_table EnergyDipoleRMSE with AtomicDipolesMACE model"
            )
        model = modules.EnergyDipolesMACE(
            **model_config,
            correlation=args.correlation,
            gate=modules.gate_dict[args.gate],
            interaction_cls_first=modules.interaction_classes[
                "RealAgnosticInteractionBlock"
            ],
            MLP_irreps=o3.Irreps(args.MLP_irreps),
        )
    else:
        raise ValueError(f"Unknown model: '{args.model}'")

    model.to(device)

    # Optimizer
    decay_interactions = {}
    no_decay_interactions = {}
    for name, param in model.interactions.named_parameters():
        if "linear.weight" in name or "skip_tp_full.weight" in name:
            decay_interactions[name] = param
        else:
            no_decay_interactions[name] = param

    param_options = dict(
        params=[
            {
                "name": "embedding",
                "params": model.node_embedding.parameters(),
                "weight_decay": 0.0,
            },
            {
                "name": "interactions_decay",
                "params": list(decay_interactions.values()),
                "weight_decay": args.weight_decay,
            },
            {
                "name": "interactions_no_decay",
                "params": list(no_decay_interactions.values()),
                "weight_decay": 0.0,
            },
            {
                "name": "products",
                "params": model.products.parameters(),
                "weight_decay": args.weight_decay,
            },
            {
                "name": "readouts",
                "params": model.readouts.parameters(),
                "weight_decay": 0.0,
            },
        ],
        lr=args.lr,
        amsgrad=args.amsgrad,
    )

    optimizer: torch.optim.Optimizer
    if args.optimizer == "adamw":
        optimizer = torch.optim.AdamW(**param_options)
    else:
        optimizer = torch.optim.Adam(**param_options)

    logger = tools.MetricsLogger(directory=args.results_dir, tag=tag + "_train")

    lr_scheduler = LRScheduler(optimizer, args)

    swa: tools.SWAContainer | None = None
    swas = [False]
    if args.swa:
        if dipole_only is not False:
            raise RuntimeError("swa for dipole fitting not implemented")
        swas.append(True)
        if args.start_swa is None:
            args.start_swa = (
                args.max_num_epochs // 4 * 3
            )  # if not set start swa at 75% of training
        log_msg = (
            f"Using stochastic weight averaging (after {args.start_swa} epochs) "
            f"with energy weight : {args.swa_energy_weight}, forces weight : "
            f"{args.swa_forces_weight}, learning rate : {args.swa_lr}"
        )
        if args.loss == "forces_only":
            raise RuntimeError("Can not select SWA with forces-only loss.")
        if args.loss == "virials":
            loss_fn_energy = modules.WeightedEnergyForcesVirialsLoss(
                energy_weight=args.swa_energy_weight,
                forces_weight=args.swa_forces_weight,
                virials_weight=args.swa_virials_weight,
            )
        elif args.loss == "stress":
            loss_fn_energy = modules.WeightedEnergyForcesStressLoss(
                energy_weight=args.swa_energy_weight,
                forces_weight=args.swa_forces_weight,
                stress_weight=args.swa_stress_weight,
            )
        elif args.loss == "energy_forces_dipole":
            loss_fn_energy = modules.WeightedEnergyForcesDipoleLoss(
                args.swa_energy_weight,
                forces_weight=args.swa_forces_weight,
                dipole_weight=args.swa_dipole_weight,
            )
            print(f"{log_msg}, dipole weight : {args.swa_dipole_weight}")
        else:
            loss_fn_energy = modules.WeightedEnergyForcesLoss(
                energy_weight=args.swa_energy_weight,
                forces_weight=args.swa_forces_weight,
            )
            print(log_msg)
        swa = tools.SWAContainer(
            model=AveragedModel(model),
            scheduler=SWALR(
                optimizer=optimizer,
                swa_lr=args.swa_lr,
                anneal_epochs=1,
                anneal_strategy="linear",
            ),
            start=args.start_swa,
            loss_fn=loss_fn_energy,
        )

    checkpoint_handler = tools.CheckpointHandler(
        directory=args.checkpoints_dir,
        tag=tag,
        keep=args.keep_checkpoints,
        swa_start=args.start_swa,
    )

    start_epoch = 0
    if args.restart_latest:
        try:
            opt_start_epoch = checkpoint_handler.load_latest(
                state=tools.CheckpointState(model, optimizer, lr_scheduler),
                swa=True,
                device=device,
            )
        except Exception:  # pylint: disable=W0703
            opt_start_epoch = checkpoint_handler.load_latest(
                state=tools.CheckpointState(model, optimizer, lr_scheduler),
                swa=False,
                device=device,
            )
        if opt_start_epoch is not None:
            start_epoch = opt_start_epoch

    ema: ExponentialMovingAverage | None = None
    if args.ema:
        ema = ExponentialMovingAverage(model.parameters(), decay=args.ema_decay)

    print(model)
    print(f"Number of parameters: {tools.count_parameters(model)}")
    print(f"{optimizer=}")

    if args.wandb:
        print("Using Weights and Biases for logging")
        import wandb

        wandb_config = {}
        args_dict = vars(args)
        args_dict_json = json.dumps(args_dict)
        for key in args.wandb_log_hypers:
            wandb_config[key] = args_dict[key]
        tools.init_wandb(
            project=args.wandb_project,
            entity=args.wandb_entity,
            name=args.wandb_name,
            config=wandb_config,
        )
        wandb.run.summary["params"] = args_dict_json

    if args.distributed:
        distributed_model = DistributedDataParallel(model, device_ids=[local_rank])
    else:
        distributed_model = None

    tools.train(
        model=model,
        loss_fn=loss_fn,
        train_loader=train_loader,
        valid_loader=valid_loader,
        optimizer=optimizer,
        lr_scheduler=lr_scheduler,
        checkpoint_handler=checkpoint_handler,
        eval_interval=args.eval_interval,
        start_epoch=start_epoch,
        max_num_epochs=args.max_num_epochs,
        logger=logger,
        patience=args.patience,
        output_args=output_args,
        device=device,
        swa=swa,
        ema=ema,
        max_grad_norm=args.clip_grad,
        log_errors=args.error_table,
        log_wandb=args.wandb,
        distributed=args.distributed,
        distributed_model=distributed_model,
        train_sampler=train_sampler,
        rank=rank,
        keep_last=True,
    )

    print("Computing metrics for training, validation, and test sets")

    all_data_loaders = {
        "train": train_loader,
        "valid": valid_loader,
    }

    test_sets = {}
    if args.train_file.endswith(".xyz"):
        for name, subset in collections.tests:
            test_sets[name] = [
                mace.data.AtomicData.from_config(
                    config, z_table=z_table, cutoff=args.r_max
                )
                for config in subset
            ]
    else:
        test_files = get_files_with_suffix(args.test_dir, "_test.h5")
        for test_file in test_files:
            name = os.path.splitext(os.path.basename(test_file))[0]
            test_sets[name] = HDF5Dataset(test_file, r_max=args.r_max, z_table=z_table)

    for test_name, test_set in test_sets.items():
        test_sampler = None
        if args.distributed:
            test_sampler = torch.utils.data.distributed.DistributedSampler(
                test_set,
                num_replicas=world_size,
                rank=rank,
                shuffle=True,
                drop_last=True,
                seed=args.seed,
            )
        test_loader = torch_geometric.dataloader.DataLoader(
            test_set,
            batch_size=args.valid_batch_size,
            shuffle=(test_sampler is None),
            drop_last=test_set.drop_last,  # type: ignore[attr-defined]
            num_workers=args.num_workers,
            pin_memory=args.pin_memory,
        )
        all_data_loaders[test_name] = test_loader

    for swa_eval in swas:
        epoch = checkpoint_handler.load_latest(
            state=tools.CheckpointState(model, optimizer, lr_scheduler),
            swa=swa_eval,
            device=device,
        )
        model.to(device)
        if args.distributed:
            distributed_model = DistributedDataParallel(model, device_ids=[local_rank])
        model_to_evaluate = model if not args.distributed else distributed_model
        print(f"Loaded model from epoch {epoch}")

        table = create_error_table(
            table_type=args.error_table,
            all_data_loaders=all_data_loaders,
            model=model_to_evaluate,
            loss_fn=loss_fn,
            output_args=output_args,
            log_wandb=args.wandb,
            device=device,
            distributed=args.distributed,
        )
        print("\n" + str(table))

        if rank == 0:
            # Save entire model
            if swa_eval:
                model_path = Path(args.checkpoints_dir) / (tag + "_swa.model")
            else:
                model_path = Path(args.checkpoints_dir) / (tag + ".model")
            print(f"Saving model to {model_path}")
            if args.save_cpu:
                model = model.to("cpu")
            torch.save(model, model_path)

            if swa_eval:
                torch.save(model, Path(args.model_dir) / (args.name + "_swa.model"))
            else:
                torch.save(model, Path(args.model_dir) / (args.name + ".model"))

        if args.distributed:
            torch.distributed.barrier()

    print("Done")
    if args.distributed:
        torch.distributed.destroy_process_group()


if __name__ == "__main__":
    main(
        name=job_name,
        train_file="../train.h5",
        valid_file="../valid.h5",
        statistics_file="../statistics.json",
        loss="stress",
        energy_weight=1,
        forces_weight=1,
        compute_stress=True,
        stress_weight=0.01,
        eval_interval=1,
        config_type_weights='{"Default":1.0}',
        E0s="average",
        error_table="PerAtomMAE",
        stress_key="stress",
        model="ScaleShiftMACE",
        MLP_irreps="64x0e",
        interaction_first="RealAgnosticResidualInteractionBlock",
        interaction="RealAgnosticResidualInteractionBlock",
        num_interactions=2,
        max_ell=3,
        hidden_irreps="64x0e + 64x1o + 64x2e",
        num_cutoff_basis=10,
        lr=0.005,
        correlation=3,
        r_max=6.0,
        num_radial_basis=10,
        scaling="rms_forces_scaling",
        distributed=True,
        num_workers=8,
        batch_size=10,
        valid_batch_size=30,
        max_num_epochs=150,
        patience=50,
        amsgrad=True,
        weight_decay=1e-8,
        ema=True,
        ema_decay=0.999,
        default_dtype="float32",
        clip_grad=100,
        device="cuda",
        seed=3,
        save_cpu=True,
        wandb=True,
        wandb_entity=WANDB_PATH.split("/")[0],
        wandb_project=WANDB_PATH.split("/")[1],
        wandb_name=job_name,
        restart_latest=True,
    )
