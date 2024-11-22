#!/bin/bash


## this script is for imitation of SLURM tasks array on local machine with multiple GPUs
export TF_CPP_MIN_LOG_LEVEL=3

THREADS=4  # number of threads for OMP, MKL, NUMEXPR etc. To share resources on single machine

export NGPU=4  # Set the total number of GPUs here

model_name="MP_GRACE_2L_r6_11Nov2024" # just for information, model name is hardcoded in 1_test_srme.py

export MODEL_NAME="${model_name}"
echo "MODEL_NAME=${MODEL_NAME}"
export MKL_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS

export SLURM_ARRAY_TASK_COUNT=96  # slurm_array_task_count

export SLURM_ARRAY_JOB_ID="production"
export SLURM_JOB_ID="production"


# Function to kill all child processes when the script receives a termination signal
cleanup() {
  echo "Terminating all child processes..."
  pkill -P $$
  exit 1
}

# Set trap to catch signals and trigger the cleanup function
trap cleanup SIGINT SIGTERM


#exec > "${MODEL_NAME}.out" 2>&1
echo "MODEL_NAME=${MODEL_NAME}"
echo "SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"


echo "Running test_grace.py"
for task_id in $(seq 0 $((SLURM_ARRAY_TASK_COUNT-1)))
do
  CUDA_VISIBLE_DEVICES=$((task_id % NGPU)) SLURM_ARRAY_TASK_ID=${task_id}  python test_grace.py > "out-${task_id}.txt" 2>&1 &
done
wait
