from __future__ import annotations

import subprocess

__author__ = "Yuan Chiang"
__date__ = "2023-09-18"


# Use multi-gpu distributed training branch
branch_repo_url = "git+https://github.com/chiang-yuan/mace.git@mbd"

try:
    subprocess.run(["module load pytorch/2.0.1"], check=True)
    subprocess.run([". ~/.venv/py311/bin/activate"], check=True)
    subprocess.run(["pip", "uninstall", "mace", "-y"], check=True)
    subprocess.run(["pip", "install", branch_repo_url], check=True)
    print("Package installed successfully.")
except Exception as e:
    print(f"Error installing the package: {e}")


job_name = "mace-train-mptrj"
account="matgen"

total_time="7-00:00:00"
time="2:00:00"
time_min="1:00:00"

sbatch_script_path="./mace-universal-distributed.sbatch"

sbatch_command = [
    "sbatch",
    "--account=" + account,
    "--job-name=" + job_name,
    "--time=" + time,
    "--time-min=" + time_min,
    "--comment=" + total_time,
    sbatch_script_path
]

try:
    subprocess.run(sbatch_command, check=True)
    print("Sbatch script executed successfully.")
except Exception as e:
    print(f"Error executing the sbatch script: {e}")
