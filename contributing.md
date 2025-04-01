# How to submit new models to Matbench Discovery

## 🔨 &thinsp; Installation

Clone the repo and install `matbench_discovery` into your Python environment (`--config-settings editable-mode=compat` helps the VS Code Python extension resolve `matbench_discovery` imports):

```zsh
git clone https://github.com/janosh/matbench-discovery --depth 1
pip install -e ./matbench-discovery --config-settings editable-mode=compat
```

There's also a [PyPI package](https://pypi.org/project/matbench-discovery) for faster installation if you don't need the latest code changes (unlikely if you're planning to submit a model since the benchmark is under active development).

## ✨ &thinsp; How to submit a new model

To submit a new model to this benchmark and add it to our leaderboard, please create a pull request to the [`main` branch][repo] that includes the 3 required items:

1. You should share your model's predictions through a cloud storage service (we recommend [Figshare](https://figshare.com)) and include the download links in your PR description. Your cloud storage directory should contain files in a compressed format with the following naming convention: `<arch-name>/<model-variant>/<yyyy-mm-dd>-<eval-task>.{csv.gz|json.gz}`. For example, a in the case of MACE-MP-0, the file paths would be:

   - geometry optimization: `mace/mace-mp-0/2023-12-11-wbm-IS2RE-FIRE.jsonl.gz` (use [JSON Lines format](https://jsonlines.org) for fast loading of small numbers of structures with `pandas.read_json(lines=True, nrows=100)` for inspection)
   - discovery: `mace/mace-mp-0/2023-12-11-wbm-IS2RE.csv.gz`
   - phonons: `mace/mace-mp-0/2024-11-09-kappa-103-FIRE-dist=0.01-fmax=1e-4-symprec=1e-5.json.gz`

   The files should contain the following information:

   1. `<arch-name>/<model-variant>/<yyyy-mm-dd>-wbm-geo-opt-<optimizer>.json.gz`: The model's relaxed structures as compressed JSON containing:

      - Final relaxed structures (as ASE `Atoms` or pymatgen `Structures`)
      - Final energies (eV), forces (eV/Å), stress (eV/Å³) and volume (Å³)
      - Material IDs matching the WBM test set

   2. `<arch-name>/<model-variant>/<yyyy-mm-dd>-wbm-IS2RE.csv.gz`: A compressed CSV file with:

      - Material IDs matching the WBM test set
      - Final formation energies per atom (eV/atom)

   3. `<arch-name>/<model-variant>/<yyyy-mm-dd>-kappa-103-FIRE-<values-of-dist|fmax|symprec>.json.gz`: A compressed JSON file with:
      - Material IDs matching the WBM test set
      - Predicted thermal conductivity (κ) values (W/mK)

   Optionally you can also share the complete relaxation trajectories. Having the forces and stresses at each relaxation step also allows analyzing any pathological behavior for structures where relaxation failed or went haywire.

   Example of how to record these quantities for a single structure with ASE:

   ```py
   from collections import defaultdict
   import pandas as pd
   from ase.atoms import Atoms
   from ase.optimize import FIRE
   from mace.calculators import mace_mp

   trajectory = defaultdict(list)
   batio3 = Atoms(
      "BaTiO3",
      scaled_positions=[
         (0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0, 0.5), (0, 0.5, 0.5), (0.5, 0.5, 0)
      ],
      cell=[4] * 3,
   )
   batio3.calc = mace_mp(model_name="medium", default_dtype="float64")

   def callback() -> None:
      """Record energy, forces, stress and volume at each step."""
      trajectory["energy"] += [batio3.get_potential_energy()]
      trajectory["forces"] += [batio3.get_forces()]
      trajectory["stress"] += [batio3.get_stress()]
      trajectory["volume"] += [batio3.get_volume()]
      # Optionally save structure at each step (results in much larger files)
      trajectory["atoms"] += [batio3.copy()]

   opt = FIRE(batio3)
   opt.attach(callback) # register callback
   opt.run(fmax=0.01, steps=500)  # optimize geometry

   df_traj = pd.DataFrame(trajectory)
   df_traj.index.name = "step"
   df_traj.to_csv("trajectory.csv.gz")  # Save final structure and trajectory data
   ```

1. `test_<model_name>_discovery.py`: The Python script that generated the WBM final energy predictions given the initial (unrelaxed) DFT structures. Ideally, this file should have comments explaining at a high level what the code is doing and how the model works so others can understand and reproduce your results. If the model deployed on this benchmark was trained specifically for this purpose (i.e. if you wrote any training/fine-tuning code while preparing your PR), please also include it as `train_<model_name>.py`.

1. `<model_name.yml>`: A file to record all relevant metadata of your algorithm like model name and version, authors, package requirements, links to publications, notes, etc. Here's a template:

   ```yml
   model_name: My new model # required (this must match the model's label which is the 3rd arg in the matbench_discovery.preds.Model enum)
   model_key: my-new-model # this should match the name of the YAML file and determines the URL /models/<model_key> on which details of the model are displayed on the website
   model_version: 1.0.0
   matbench_discovery_version: 1.0
   date_added: "2023-01-01"
   date_published: "2022-12-05"
   authors:
     - name: John Doe
       affiliation: Some University, Some National Lab
       email: john-doe@uni.edu
       orcid: https://orcid.org/0000-xxxx-yyyy-zzzz
       url: lab.gov/john-doe
       corresponding: true
       role: Model & PR
     - name: Jane Doe
       affiliation: Some National Lab
       email: jane-doe@lab.gov
       url: uni.edu/jane-doe
       orcid: https://orcid.org/0000-xxxx-yyyy-zzzz
       role: Model
   repo: https://github.com/<user>/<repo>
   url: https://<model-docs-or-similar>.org
   doi: https://doi.org/10.5281/zenodo.0000000
   preprint: https://arxiv.org/abs/xxxx.xxxxx
   checkpoint_url: https://github.com/<org>/<repo>/releases/download/<tag>/<model-file> # checkpoint URL can be any URL but GitHub release URLs are encouraged

   openness: OSOD # see `Open` enum in matbench_discovery/enums.py
   train_task: S2EFS # see `Task` enum in matbench_discovery/enums.py
   test_task: IS2RE-SR # see `Task` enum in matbench_discovery/enums.py
   targets: EFS_G # see `Targets` enum in matbench_discovery/enums.py
   model_type: UIP # see `ModelType` enum in matbench_discovery/enums.py
   model_params: 5_000_000
   trained_for_benchmark: true
   n_estimators: 1

   hyperparams: # strongly recommended to list relaxation hyperparams
      max_force: 0.05
      max_steps: 500
      ase_optimizer: FIRE
      optimizer: Adam
      ... # additional hyperparameters describing training

   training_cost: # list any hardware used to train the model and for how long
     # hardware: { amount: float, hours: float, cost: float [USD] }
     Nvidia A100 GPUs: { amount: 4, hours: 30, cost: 1000 }
     Nvidia H100 GPUs: { amount: 1, hours: 90, cost: 2000 }
     Intel Xeon Gold 6254 CPUs: { amount: 12, hours: 100, cost: 100 }

   requirements: # strongly recommended
     torch: 1.13.0
     torch-geometric: 2.0.9
     ...

   training_set: [MPtrj] # list of keys from data/training-sets.yml

   notes: # notes can have any key, be multiline and support markdown.
     description: This is how my model works...
     steps: |
      Optional *free-form* [markdown](example.com) notes.

   metrics:
     phonons:
       kappa_103:
         pred_file: models/<model_dir>/<yyyy-mm-dd>-kappa-103-<values-of-dist|fmax|symprec>.json.gz
         pred_file_url: https://ndownloader.figshare.com/files/<figshare_id>
     geo_opt: # only applicable if the model performed structure relaxation
       pred_file: models/<model_dir>/<yyyy-mm-dd>-wbm-geo-opt-<optimizer>.json.gz # should contain the models relaxed structures as ASE Atoms or pymatgen Structures, and separate columns for material_id and energies/forces/stresses at each relaxation step
       pred_file_url: https://ndownloader.figshare.com/files/<figshare_id>
       struct_col: <column_name_of_material_ids_in_relaxed_structures>
     discovery:
       pred_file: models/<model_dir>/<yyyy-mm-dd>-<model_name>-wbm-IS2RE.csv.gz # should contain the models energy predictions for the WBM test set
       pred_file_url: https://ndownloader.figshare.com/files/<figshare_id>
       pred_col: e_form_per_atom_<model_name>
   ```

   Arbitrary other keys can be added as needed. The above keys will be schema-validated with `pre-commit` (if installed) with errors for missing keys.

Please see any of the subdirectories in [`models/`](https://github.com/janosh/matbench-discovery/bob/-/models) for example submissions. More detailed step-by-step instructions below.

### Step 1: Clone the repo

```sh
git clone https://github.com/janosh/matbench-discovery --depth 1
cd matbench-discovery
git checkout -b model-name-you-want-to-add
```

Tip: `--depth 1` only clones the latest commit, not the full `git history` which is faster if a repo contains large data files that changed over time.

### Step 2: Commit model metadata and train/test scripts to the repo

Create a new folder

```sh
mkdir models/<model_name>
```

and place the above-listed files there. The file structure should look like this:

```txt
matbench-discovery-root
└── models
    └── <model_name>
        ├── <model_name>.yml
        ├── test_<model_name>_discovery.py
        ├── test_<model_name>_kappa.py
        ├── readme.md  # optional
        └── train_<model_name>.py  # optional
```

You can include arbitrary other supporting files like metadata and model features (below 10MB total to keep `git clone` time low) if they are needed to run the model or help others reproduce your results. For larger files, please upload to [Figshare](https://figshare.com) or similar and share the link in your PR description.

Add the model to the `Model` enum in [`matbench_discovery/enums.py`](https://github.com/janosh/matbench-discovery/blob/57d0d0c8a14cd3/matbench_discovery/enums.py#L274) pointing to the correct metadata file.

### Step 3: Upload results files to Zenodo

Upload the files for the predictions, optimized geometries and phonons to Zenodo and share the links in your PR description.

### Step 4: Open a PR to the [Matbench Discovery repo][repo]

Commit your files to the repo on a branch called `<model_name>` and create a pull request (PR) to the Matbench repository.

```sh
git add -a models/<model_name>
git commit -m 'add <model_name> to Matbench Discovery leaderboard'
```

And you're done! Once tests pass and the PR is merged, your model will be added to the leaderboard! 🎉

### Step 5 (optional): Copy WandB runs into our project

[Weights and Biases](https://wandb.ai) is a tool for logging training and test runs of ML models. It's free, (partly) [open source](https://github.com/wandb/wandb) and offers a [special plan for academics](https://wandb.ai/site/research). It auto-collects metadata like

- what hardware the model is running on
- and for how long,
- what the CPU, GPU and network utilization was over that period,
- the exact code in the script that launched the run, and
- which versions of dependencies were installed in the environment your model ran in.

This information can be useful for others looking to reproduce your results or compare their model to yours i.t.o. computational cost. We therefore strongly recommend tracking all runs that went into a model submission with WandB so that the runs can be copied over to our WandB project at <https://wandb.ai/janosh/matbench-discovery> for everyone to inspect. This also allows us to include your model in more detailed analysis (see the SI in the [preprint](https://arxiv.org/abs/2308.14920)).

## 😵‍💫 &thinsp; Troubleshooting

Having problems? Please [open an issue on GitHub](https://github.com/janosh/matbench-discovery/issues). We're happy to help! 😊

[repo]: https://github.com/janosh/matbench-discovery
