#!/bin/bash
#SBATCH -q preempt
#SBATCH -C gpu
#SBATCH -G 40
#SBATCH -N 10
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
##SBATCH --time-min=01:00:00
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
##SBATCH --comment=7-00:00:00
#SBATCH --signal=B:USR1@180
#SBATCH --requeue
#SBATCH --exclusive
#SBATCH --open-mode=append


module load pytorch/2.0.1
. ~/.venv/py311/bin/activate

exp_name=$(basename "$SLURM_SUBMIT_DIR")

srun python ./run_train.py \
    --name=$exp_name \
    --train_file="../train.h5" \
    --valid_file="../valid.h5" \
    --statistics_file="../statistics.json" \
    --loss='uip' \
    --energy_weight=1 \
    --forces_weight=1 \
    --compute_stress=True \
    --stress_weight=0.01 \
    --eval_interval=1 \
    --config_type_weights='{"Default":1.0}' \
    --E0s='average' \
    --error_table='PerAtomMAE' \
    --stress_key='stress' \
    --model="ScaleShiftMACE" \
    --MLP_irreps="64x0e" \
    --interaction_first="RealAgnosticResidualInteractionBlock" \
    --interaction="RealAgnosticResidualInteractionBlock" \
    --num_interactions=2 \
    --num_channels=128 \
    --max_ell=3 \
    --hidden_irreps='64x0e + 64x1o + 64x2e' \
    --num_cutoff_basis=10 \
    --lr=0.005 \
    --correlation=3 \
    --r_max=6.0 \
    --num_radial_basis=10 \
    --scaling='rms_forces_scaling' \
    --distributed \
    --num_workers=8 \
    --batch_size=10 \
    --valid_batch_size=30 \
    --max_num_epochs=150 \
    --patience=50 \
    --amsgrad \
    --weight_decay=1e-8 \
    --ema \
    --ema_decay=0.999 \
    --default_dtype="float32"\
    --clip_grad=100 \
    --device=cuda \
    --seed=3 \
    --save_cpu \
    --wandb --wandb_project <project> --wandb_entity <entity> --wandb_name $exp_name \
    --restart_latest &

ckpt_command=

. ./setup.sh
requeue_job func_trap USR1

wait
