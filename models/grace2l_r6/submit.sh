#!/bin/bash
#SBATCH --job-name=matbench-test
#SBATCH --time=1-00:00:00          # Maximum runtime (1 day)
#SBATCH --partition=shorttime      # Partition name
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=7                 # One task per job
#SBATCH --cpus-per-task=4          # 4 CPUs per job
#SBATCH --array=1-4               # Job array with tasks (adjust as needed)
#SBATCH --output=output_%A_%a.log  # Output file for each job
#SBATCH --exclusive                # Ensure only these jobs run on the node
#SBATCH --mem-per-cpu=6G            # Memory per CPU (adjust as needed)


# Print some job details
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "SLURM Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Using ${SLURM_CPUS_PER_TASK} CPUs on node(s): ${SLURM_NODELIST}"

uss=$(whoami)
find /dev/shm/ -user $uss -type f -mmin +30 -delete
source $HOME/.bashrc
micromamba activate grace


THREADS=${SLURM_CPUS_PER_TASK}  # number of threads for OMP, MKL, NUMEXPR etc. To share resources on single machine
export MKL_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS

export TF_NUM_INTEROP_THREADS=$THREADS
export TF_NUM_INTRAOP_THREADS=$THREADS

python test_grace.py

