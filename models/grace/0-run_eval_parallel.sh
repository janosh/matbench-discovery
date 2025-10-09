#!/bin/bash

########## HOW-TO ############
#
# STAGE 1
# 1. set GPUS to list of available GPUs
# 2. set model_name=....
# 3. adjust SLURM_ARRAY_TASK_COUNT= if needed
# 4. run locally or submit to the queue
#
# STAGE 2
# python join_grace_preds.py  2025-02-06-GRACE-1L-OAM_2Feb25-wbm-IS2RE-FIRE  # or other path
#
# STAGE 3
# upload large files to fileshare, commit models' YAML files
################################



## this script is for imitation of SLURM tasks array on local machine with multiple GPUs
export TF_CPP_MIN_LOG_LEVEL=3

THREADS=4  # number of threads for OMP, MKL, NUMEXPR etc. To share resources on single machine

# Define the list of available GPUs
GPUS=(0 1 2 3)  # Set the specific GPU IDs you want to use here
NGPU=${#GPUS[@]} # Calculate the number of GPUs based on the array length

#model_name="MP_GRACE_2L_r6_11Nov2024"
#model_name="GRACE-1L-OAM_2Feb25"
model_name="GRACE-2L-OAM-L"

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
  # Calculate the GPU index using the modulo operator (%)
  gpu_index=$((task_id % NGPU))
  # Get the actual GPU ID from the GPUS array
  gpu_id=${GPUS[$gpu_index]}

  CUDA_VISIBLE_DEVICES=$gpu_id SLURM_ARRAY_TASK_ID=${task_id}  python 0-test_grace_discovery_new.py "${MODEL_NAME}" > "out-${task_id}.txt" 2>&1 &
done
wait
