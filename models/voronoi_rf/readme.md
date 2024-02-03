# Voronoi Random Forest

## Model Architecture

Voronoi tessellation with `matminer` featurization piped into `scikit-learn` [`RandomForestRegressor`](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor).

## Reasoning

The idea behind this combination of features and model was to have an easy-to-implement baseline algorithm. It's a bit outdated in that it uses handcrafted Magpie features (which have been shown to underperform learned features on datasets exceeding ~10^4 samples) but not so weak as to be indefensible. The fact that its featurization based on Voronoi tessellation is invariant to crystal structure relaxation makes it a natural choice for predicting unrelaxed crystal stability.

## OOM errors during featurization

There was an obstacle that actually made this model more difficult to train and test than anticipated. `matminer` uses `multiprocessing` which seems to be the cause of out-of-memory errors on large structures. Initially couldn't get [`MultipleFeaturizer`] to run without crashing even when running on small subsets of the data (1%) and setting `sbatch` flag `--mem 100G`:

```log
MultipleFeaturizer:  28%|██▊       | 724/2575 [01:08<04:15,  7.25it/s]/var/spool/slurm/slurmd/job7401930/slurm_script: line 4: 2625851 Killed                  python
slurmstepd: error: Detected 52 oom-kill event(s) in StepId=7401930.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
4:00
```

Saving tip came from [Alex Dunn via Slack](https://berkeleytheory.slack.com/archives/D03ULSTNRMX/p1668746161675349) to try `featurizer.set_n_jobs(1)`.

## Archive

Files in `2022-10-04-rhys-voronoi.zip` received from Rhys via [Slack](https://ml-physics.slack.com/archives/DD8GBBRLN/p1664929946687049). They are unchanged originals.

[`multiplefeaturizer`]: https://hackingmaterials.lbl.gov/matminer/matminer.featurizers.html#matminer.featurizers.base.MultipleFeaturizer
