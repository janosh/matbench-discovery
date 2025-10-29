#!/bin/bash

########## HOW-TO ############
#
# Stage 1
#   1. Set GPUS to list of available GPUs
#   2. Set model_name variable below
#   3. Adjust SLURM_ARRAY_TASK_COUNT if needed
#   4. Run locally or submit to queue
#
# Stage 2
#   python 2-join_matbench_preds.py
#
# Stage 3
#   python 3-filter_bad_preds.py
#
# Stage 4
#   Upload large files to figshare, commit models' YAML files
################################

## This script imitates SLURM task arrays on a local machine with multiple GPUs
export TF_CPP_MIN_LOG_LEVEL=3

THREADS=4  # Number of threads for OMP, MKL, NUMEXPR etc. to share resources on single machine

# Define the list of available GPUs
GPUS=(0 1 2 3)  # Set the specific GPU IDs you want to use here
NGPU=${#GPUS[@]}  # Calculate the number of GPUs based on the array length

export MODEL_NAME="GRACE-2L-OAM-L"
export MKL_NUM_THREADS=${THREADS}
export NUMEXPR_NUM_THREADS=${THREADS}
export OMP_NUM_THREADS=${THREADS}

export SLURM_ARRAY_TASK_COUNT=96

export SLURM_ARRAY_JOB_ID="production"
export SLURM_JOB_ID="production"


# Function to kill all child processes when the script receives a termination signal
cleanup() {
  echo "Terminating all child processes..."
  kill -- -$$ 2>/dev/null  # Kill the entire process group
  exit 1
}

# Set trap to catch signals and trigger the cleanup function
trap cleanup SIGINT SIGTERM

echo "MODEL_NAME=${MODEL_NAME}"
echo "SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "Running 1-test_grace_discovery.py"

for task_id in $(seq 0 $((SLURM_ARRAY_TASK_COUNT - 1))); do
  # Calculate the GPU index using the modulo operator
  gpu_index=$((task_id % NGPU))
  # Get the actual GPU ID from the GPUS array
  gpu_id=${GPUS[$gpu_index]}

  CUDA_VISIBLE_DEVICES=${gpu_id} SLURM_ARRAY_TASK_ID=${task_id} python 1-test_grace_discovery.py > "out-${task_id}.txt" 2>&1 &
done
wait
