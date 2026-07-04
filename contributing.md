# How to submit new models to Matbench Discovery

## 🔨 &thinsp; Installation

Clone [the repo](https://github.com/janosh/matbench-discovery) and install `matbench-discovery` into your Python environment:

```zsh
git clone https://github.com/janosh/matbench-discovery --depth 1
pip install -e ./matbench-discovery
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

   ### Script Dependencies Declaration (Required)

   **All test scripts must include a script dependencies section** using the [PEP 723 inline metadata format](https://docs.astral.sh/uv/guides/scripts/#declaring-script-dependencies). This is critical for reproducibility and allows others to easily rerun your code with `uv run test_<model_name>_<task>.py` without manually setting up virtual environments or installing dependencies.

   #### Why This Matters

   - **Reproducibility**: Future readers can still run your scripts even if there are breaking changes in package releases
   - **Convenience**: One command (`uv run script.py`) installs dependencies and runs the script
   - **Documentation**: Makes it clear which exact package versions were used for your submission

   ```python
   # /// script
   # requires-python = ">=3.11,<3.13"
   # dependencies = [
   #   "numpy>=2.3.4",
   #   "torch>=2.8.0",
   #   "ase>=3.26.0",
   #   "pymatgen>=2025.10.7",
   #   "matbench-discovery>=1.3.1",
   # ]
   #
   # [tool.uv.sources]
   # matbench-discovery = { path = "../../", editable = true }
   # ///

   import numpy as np
   import torch
   from ase.atoms import Atoms
   # ... rest of your script
   ```

   ### Automatic Checkpoint Downloading (Recommended)

   Your script should automatically download the model checkpoint when first run and cache it to a standard location such as `~/.cache/<model_name>/checkpoint.pt` for later reuse. This also lowers the barrier for others to reproduce results.

1. `<model_name.yml>`: A file to record all relevant metadata of your algorithm like model name and version, authors, package requirements, links to publications, notes, etc. Here's a template:

   ```yml
   model_name: My new model # required (this must match the model's label which is the 3rd arg in the matbench_discovery.preds.Model enum)
   model_key: my-new-model # this should match the name of the YAML file and determines the URL /models/<model_key> on which details of the model are displayed on the website
   model_version: 1.0.0
   date_added: '2023-01-01'
   date_published: '2022-12-05'
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
   checkpoint_url: https://url.of/model-checkpoint-file # should be the direct download URL of the model checkpoint file. If there's a page documenting the checkpoint, append that as a YAML comment.

   openness: OSOD # see `Open` enum in matbench_discovery/enums.py
   train_task: S2EFS # see `Task` enum in matbench_discovery/enums.py
   test_task: IS2RE-SR # see `Task` enum in matbench_discovery/enums.py
   targets: EFS_G # see `Targets` enum in matbench_discovery/enums.py
   model_type: UIP # see `ModelType` enum in matbench_discovery/enums.py
   model_params: 5_000_000
   trained_for_benchmark: true
   n_estimators: 1

   license:
     code: MIT
     code_url: https://url.of/code-license
     checkpoint: CC-BY-4.0
     # URL that points to the license file for the model checkpoint, not the checkpoint file itself.
     checkpoint_url: https://url.of/model-checkpoint-license

   hyperparams: # strongly recommended to list relaxation hyperparams
     max_force: 0.05
     max_steps: 500
     ase_optimizer: FIRE
     optimizer: Adam
     graph_construction_radius: 6.0
     max_neighbors: 50
       ... # additional hyperparameters describing training and inference

   training_cost: # list any hardware used to train the model and for how long
     # hardware: { amount: float, hours: float, cost: float [USD] }
     Nvidia A100 GPUs: { amount: 4, hours: 30, cost: 1000 }
     Nvidia H100 GPUs: { amount: 1, hours: 90, cost: 2000 }
     Intel Xeon Gold 6254 CPUs: { amount: 12, hours: 100, cost: 100 }

   requirements: # strongly recommended
     torch: 1.13.0
     torch-geometric: 2.0.9
       ...

   training_set: [MPtrj] # list of keys from data/datasets.yml

   notes: # notes can have any key, be multiline and support markdown.
     description: This is how my model works...
     steps: |
       Optional *free-form* [markdown](example.com) notes.

   metrics:
     phonons:
       kappa_103:
         pred_file: models/<model_dir>/<yyyy-mm-dd>-kappa-103-<values-of-dist|fmax|symprec>.json.gz
         pred_file_url: https://figshare.com/files/<figshare_id>
     geo_opt: # only applicable if the model performed structure relaxation
       pred_file: models/<model_dir>/<yyyy-mm-dd>-wbm-geo-opt-<optimizer>.json.gz # should contain the models relaxed structures as ASE Atoms or pymatgen Structures, and separate columns for material_id and energies/forces/stresses at each relaxation step
       pred_file_url: https://figshare.com/files/<figshare_id>
       struct_col: <column_name_of_material_ids_in_relaxed_structures>
     discovery:
       pred_file: models/<model_dir>/<yyyy-mm-dd>-<model_name>-wbm-IS2RE.csv.gz # should contain the models energy predictions for the WBM test set
       pred_file_url: https://figshare.com/files/<figshare_id>
       pred_col: e_form_per_atom_<model_name>
     md: # optional: finite-temperature molecular dynamics metrics
       pred_file: models/<model_dir>/<yyyy-mm-dd>-<model_name>-md-metrics.csv.gz # per-system energy/force RMSE, RDF/vDOS errors and pressure metrics vs ab-initio reference trajectories
       pred_file_url: https://figshare.com/files/<figshare_id>
   ```

   Arbitrary other keys can be added as needed. The above keys will be schema-validated with `prek` (if installed) with errors for missing keys.

Please see any of the subdirectories in [`models/`](https://github.com/janosh/matbench-discovery/blob/-/models) for example submissions. More detailed step-by-step instructions below.

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

The molecular dynamics task does not need a per-model script. Instead register your model's ASE calculator and its `uv` dependencies once in the shared calculator registry [`matbench_discovery/calculators.py`](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/calculators.py) (the same registry also powers the diatomics task via `models/run_diatomics.py`), then run it through the shared runner (which auto-downloads the label-free DynaMat v1.0 reference set):

```sh
# smoke-test the pipeline in seconds, then launch the full 20 ps NVT benchmark
uv run --with <your-deps> models/run_md.py --model <model_key> --dry-run
uv run --with <your-deps> models/run_md.py --model <model_key> --write-yaml
# or let the runner print the exact uv command for your model:
uv run models/run_md.py --print-cmd --model <model_key>
```

`--write-yaml` records the model-level MD metrics under `metrics.md` in your model's YAML and writes the per-system predictions to a gzipped CSV named `<yyyy-mm-dd>-<model_name>-md-metrics.csv.gz`. Public runs compute the observable metrics only; energy/force RMSEs are maintainer-computed private-label diagnostics and are not required in external submissions. Upload that CSV to Figshare (or similar) and set its download URL as the `pred_file_url` field of `metrics.md` — see [Step 3](#step-3-upload-results-files-to-figshare-or-similar) for the upload conventions and YAML field definitions.

Add the model to the `Model` enum in [`matbench_discovery/enums.py`](https://github.com/janosh/matbench-discovery/blob/57d0d0c8a14cd3/matbench_discovery/enums.py#L274) pointing to the correct metadata file.

> [!WARNING]
> Do not include a `filter_bad_preds.py` script or equivalent. This was a historical flaw in the evaluation and now all OOM outlier prediction handling is handled internally and consistently removing the need for this filtering stage.

### Step 3: Upload results files to Figshare (or similar)

Upload the files for the predictions, optimized geometries and phonons to [Figshare](https://figshare.com) (or any cloud storage with stable direct-download links) and record the download URLs in the `pred_file_url` keys of your model's YAML file. These URLs are what our automated ingestion (see below) downloads your predictions from, so they must resolve without authentication. No need to worry about link rot long-term: during ingestion we archive copies of your files to the project's own Figshare articles and update the URLs to point at them.

### Step 4: Open a PR to the [Matbench Discovery repo][repo]

Commit your files to the repo on a branch called `<model_name>` and create a pull request (PR) to the Matbench repository.

```sh
git add -a models/<model_name>
git commit -m 'add <model_name> to Matbench Discovery leaderboard'
```

And you're done! Once tests pass, a maintainer triggers automated ingestion on your PR (see below) and after merge your model appears on the leaderboard and in all site figures! 🎉

> [!IMPORTANT]
> Keep the **"Allow edits by maintainers"** checkbox on your PR enabled (it's on by default). The ingestion bot pushes evaluation results and updated site assets as a commit onto your PR branch — without that permission, ingestion has to be run manually by a maintainer.

### Step 5: Validate your submission with the justfile

We provide a [`justfile`](https://github.com/casey/just) at the repository root to help validate your submission before opening a PR. The `prepare-model-submission` command runs all evaluation scripts, generates required figures, and checks the PR checklist requirements.

#### Installing just

If you don't have `just` installed, you can install it with:

```sh
# macOS
brew install just

# Linux (via cargo)
cargo install just

# Or with pipx
pipx install rust-just

# Or see https://github.com/casey/just#installation for more options
```

#### Running the validation

```sh
just prepare-model-submission <model_name>
```

For example:

```sh
just prepare-model-submission mace-mpa-0
# or with underscores (both work)
just prepare-model-submission mace_mpa_0
```

This command will:

1. **Check PR requirements**: Verify your YAML file exists, model is registered in the `Model` enum, prediction URLs are in the YAML, and test scripts exist
2. **Run evaluation scripts**: Execute discovery, kappa (phonon), and diatomic metrics evaluation for your model
3. **Generate figures**: Create the required energy parity and per-element error plots

### Step 6 (optional): Copy WandB runs into our project

[Weights and Biases](https://wandb.ai) ([GitHub](https://github.com/wandb/wandb)) is a tool for logging training and test runs of ML models. It auto-collects metadata like:

- what hardware the model is running on
- and for how long,
- what the CPU, GPU and network utilization was over that period,
- the exact code in the script that launched the run, and
- which versions of dependencies were installed in the environment your model ran in.

This information can be useful for others looking to reproduce your results or compare their model to yours i.t.o. computational cost. We therefore strongly recommend tracking all runs that went into a model submission with WandB so that the runs can be copied over to our WandB project at <https://wandb.ai/janosh/matbench-discovery> for everyone to inspect.

## 🤖 &thinsp; Automated ingestion (what happens after your PR is opened)

Once your submission PR looks ready, a maintainer applies the `ingest-model` label. CI ([`update-site-figs.yml`](https://github.com/janosh/matbench-discovery/blob/-/.github/workflows/update-site-figs.yml)) then runs the full pipeline — no manual work needed from you:

1. **Safe checkout** — CI checks out trusted `main` code and overlays *only your PR's data*: `models/**` plus your `Model` enum addition, which is accepted only if the diff is provably additive member lines (anything else fails closed and a maintainer ingests locally after review). Submission code is never executed in CI, which is what makes it safe to run this on unreviewed PRs with credentials present.
2. **Evals + checklist** — your prediction files are downloaded from the `pred_file_url` links in your YAML, all evals run (discovery, kappa, geo-opt, diatomics), metrics are written into your model YAML, and the PR checklist is enforced.
3. **Figshare archival** — your prediction + analysis files are re-uploaded to the project's own Figshare articles (one per prediction task) for longevity, and your YAML's `*_url` keys are rewritten to the archived copies.
4. **Per-model site assets** — energy/kappa parity assets (published to the GitHub release the site build downloads from), parity manifests and per-element error data are generated.
5. **Multi-model figures** — all site figure payloads (`site/src/figs/*.json.gz`) are regenerated so every page (`/models/tmi`, `/tasks/geo-opt`, data pages) includes your model, validated by payload shape tests before committing.
6. **One commit onto your PR** — the updated YAML, payloads and manifests are pushed to your PR branch (this is why "Allow edits by maintainers" must stay enabled), retriggering CI on the result. Merging the PR lands your model fully integrated.

A post-merge + weekly fallback job regenerates the figure payloads and opens a follow-up PR only if they're stale, so the site can't silently fall out of date even if a PR merges without the label.

Maintainer notes: ingestion requires the repo secrets `SITE_FIGS_PAT` (classic PAT with `public_repo` scope, used to push to fork branches) and `FIGSHARE_TOKEN` (archival uploads). Re-trigger by re-applying the label, or run `just ingest-model <model_name>` locally for submissions whose enum diff fails validation. `just update-site-figs` refreshes the multi-model figure payloads on their own.

If CI flags the figure payloads as stale for your PR (e.g. `test_discovery_payload_covers_active_models` fails after adding or superseding a model), run `just update-site-figs <model_name>` and commit the changed `site/src/figs/*.json.gz` files. This splices only your model's freshly computed entries into the committed payloads, pulling your prediction files from the `pred_file_url` entries in the model YAML, the same links used by automated ingestion.

## 😵‍💫 &thinsp; Troubleshooting

Having problems? Please [open an issue on GitHub](https://github.com/janosh/matbench-discovery/issues). We're happy to help! 😊

[repo]: https://github.com/janosh/matbench-discovery
