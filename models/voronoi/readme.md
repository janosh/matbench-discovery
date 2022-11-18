# Voronoi Tessellation with matminer featurezation piped into `scikit-learn` Random Forest

## OOM errors during featurization

`multiprocessing` seems to be the cause of out-of-memory errors on large structures. Initially couldn't get the `matminer` `MultipleFeaturizer` to run without crashing even when running on small subsets of the data (1%) and setting `sbatch` flag `--mem 100G`:

```log
MultipleFeaturizer:  28%|██▊       | 724/2575 [01:08<04:15,  7.25it/s]/var/spool/slurm/slurmd/job7401930/slurm_script: line 4: 2625851 Killed                  python
slurmstepd: error: Detected 52 oom-kill event(s) in StepId=7401930.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
4:00
```

Saving tip came from [Alex Dunn via Slack](https://berkeleytheory.slack.com/archives/D03ULSTNRMX/p1668746161675349) to set `featurizer.set_n_jobs(1)`.

## Archive

Files in `2022-10-04-rhys-voronoi.zip` received from Rhys via [Slack](https://ml-physics.slack.com/archives/DD8GBBRLN/p1664929946687049). All originals before making any changes for this project.
