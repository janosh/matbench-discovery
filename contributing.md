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

1. You should share your model's predictions through a cloud storage service (we recommend [Figshare](https://figshare.com)) and include the stable direct-download links in your model YAML. Each artifact lives directly under `models/<family>/<model_key>/` and follows the canonical filename grammar:

   - geometry optimization: `<yyyy-mm-dd>-geo-opt.jsonl.gz` — final relaxed structures and WBM material IDs in [JSON Lines format](https://jsonlines.org). Optional analysis files use `<yyyy-mm-dd>-geo-opt-symprec=<symprec>-moyo=<version>.csv.gz`.
   - discovery: `<yyyy-mm-dd>-discovery.csv.gz` — compressed CSV with `material_id` and final `e_form_per_atom` columns.
   - phonons: `<yyyy-mm-dd>-phonons-kappa-103.json.gz` — predictions for exactly the 103 PhononDB material IDs. Optional force sets and run provenance use `-forces.json.gz` and `-run-info.json`.
   - molecular dynamics: `<yyyy-mm-dd>-md-metrics.csv.gz`.
   - diatomics: `<yyyy-mm-dd>-diatomics.json.gz`.

   Paths are validated by `matbench_discovery.data.parse_artifact_filename`. Record each
   artifact as a nested object:

   ```yaml
   pred_file:
     name: models/<family>/<model_key>/<yyyy-mm-dd>-discovery.csv.gz
     url: https://figshare.com/files/<id>
     size: 1234567  # optional
     md5: <32-char hex>  # optional; required together with size
   ```

   The same shape is used for `analysis_file`, `force_file`, and `run_info_file`.

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

2. Reproducible inference code. If your model exposes an ASE calculator, register its constructor in [`matbench_discovery/calculators.py`](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/calculators.py) and declare its isolated `uv` dependencies in YAML `environment`. The shared runners [`models/run_discovery.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_discovery.py), [`models/run_kappa.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_kappa.py), [`models/run_md.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_md.py), and [`models/run_diatomics.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_diatomics.py) are the only entry points for calculator-backed models; do not add per-model forks or join scripts.

   Calculator integration belongs in the same model PR. Normal secretless PR CI validates that code, while automated ingestion evaluates the submitted artifacts using trusted base-branch code and never executes PR code with credentials.

   Models without an ASE calculator may need a task-specific `test_<model_key>_<task>.py` or an archived direct-prediction submission — discuss this before opening the PR. If the model was trained specifically for this benchmark, include `train_<model_key>.py` or a family-level `train_<family>.py`.

   ### Script Dependencies Declaration (Required)

   **All custom test scripts must include a script dependencies section** using the [PEP 723 inline metadata format](https://docs.astral.sh/uv/guides/scripts/#declaring-script-dependencies): a single `uv run test_<model_key>_<task>.py` then installs dependencies and runs the script without manual environment setup, documents the exact package versions used for your submission, and keeps your code runnable even after breaking changes in package releases. If the resolver cannot express a required installation step such as `--no-deps`, include a documented install script as well. Standard calculator-backed models declare dependencies once in `environment` instead.

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

   Your script should automatically download the model checkpoint when first run and cache it to a standard location such as `~/.cache/<model_key>/checkpoint.pt` for later reuse. This also lowers the barrier for others to reproduce results.

3. `models/<family>/<model_key>.yml`: A file recording the model's canonical identity, release, provenance, environments, task coverage, and artifacts. It is validated by [`tests/model-schema.yml`](https://github.com/janosh/matbench-discovery/blob/main/tests/model-schema.yml). Copy a current model such as `models/mace/mace-mp-0.yml`; this abbreviated template highlights the required structure:

   ```yml
   model_name: My new model # human-readable label
   model_key: my-new-model # lowercase kebab-case URL slug (dots allowed in versions)
   model_version: 1.0.0 # upstream package / checkpoint version
   # Family is the parent models/<family>/ directory (lowercase snake_case)
   dates:
     benchmark_added: '2023-01-01'
     paper_published: '2022-12-05'
   lifecycle: active # active, superseded, deprecated, or aborted
   authors:
     - name: John Doe
       affiliation: Some University, Some National Lab
       email: john-doe@uni.edu
       orcid: https://orcid.org/0000-xxxx-yyyy-zzzz
       url: lab.gov/john-doe
       corresponding: true
     - name: Jane Doe
       affiliation: Some National Lab
       email: jane-doe@lab.gov
       url: uni.edu/jane-doe
       orcid: https://orcid.org/0000-xxxx-yyyy-zzzz
   repo: https://github.com/<user>/<repo>
   docs: https://<model-docs-or-similar>.org
   doi: https://doi.org/10.5281/zenodo.0000000
   paper: https://arxiv.org/abs/xxxx.xxxxx
   pypi: https://pypi.org/project/<package>
   pr_url: https://github.com/janosh/matbench-discovery/pull/<number>
   checkpoint_url: https://url.of/model-checkpoint-file # should be the direct download URL of the model checkpoint file. If there's a page documenting the checkpoint, append that as a YAML comment.

   license:
     code: MIT
     code_url: https://url.of/code-license
     checkpoint: CC-BY-4.0
     # URL that points to the license file for the model checkpoint, not the checkpoint file itself.
     checkpoint_url: https://url.of/model-checkpoint-license

   environment: # provenance pins + (for calculators) isolated runner install
     dependencies:
       - my-model-package>=1.0
       - torch==2.8.0
       - ase==3.26.0
       - pymatgen==2025.10.7
     python_version: '3.12'
     find_links: [] # optional
     extra_index_urls: [] # optional
     project: null # optional

   openness: OSOD
   architecture_types: [gnn]
   train_task: S2EFS
   test_task: IS2RE-SR
   targets: EFS_G
   training_sets: [MPtrj] # keys from data/datasets.yml
   # size: 10M  # optional
   # fine_tune: true  # optional
   model_params: 5_000_000
   # n_estimators: 10  # optional; omit for single models (default 1)

   # Optional — omit entirely if training hardware/cost is unknown.
   training_cost:
     entries:
       - hardware: NVIDIA A100
         count: 4
         hours_per_device: 30

   hyperparams:
     evaluation:
       max_force: 0.05
       max_steps: 500
       ase_optimizer: FIRE
       cell_filter: FrechetCellFilter
       kappa:
         protocol: phonondb-v1
         save_forces: true
     architecture:
       graph_construction_radius: 6.0
       max_neighbors: 50
     training:
       optimizer: Adam

   notes: # notes can have any key, be multiline and support markdown.
     description: This is how my model works...
     steps: |
       Optional *free-form* [markdown](example.com) notes.

   metrics:
     phonons:
       kappa_103:
         pred_file:
           name: models/<family>/<model_key>/<yyyy-mm-dd>-phonons-kappa-103.json.gz
           url: https://figshare.com/files/<figshare_id>
     geo_opt: # only applicable if the model performed structure relaxation
       pred_file:
         name: models/<family>/<model_key>/<yyyy-mm-dd>-geo-opt.jsonl.gz
         url: https://figshare.com/files/<figshare_id>
     discovery:
       pred_file:
         name: models/<family>/<model_key>/<yyyy-mm-dd>-discovery.csv.gz
         url: https://figshare.com/files/<figshare_id>
     md:
       pred_file:
         name: models/<family>/<model_key>/<yyyy-mm-dd>-md-metrics.csv.gz
         url: https://figshare.com/files/<figshare_id>
     diatomics:
       pred_file:
         name: models/<family>/<model_key>/<yyyy-mm-dd>-diatomics.json.gz
         url: https://figshare.com/files/<figshare_id>
   ```

   Task coverage is derived (no top-level `tasks:` block): metrics with results →
   `complete`; `targets: E` → force tasks `not_applicable`; `lifecycle: aborted` →
   `not_available`; otherwise `pending`. Optional exception notes live under the
   metrics task, e.g. `metrics.md: {status: pending, reason: "…"}`.

   The schema rejects unknown keys; put free-form context under `notes`.

Please see any of the subdirectories in [`models/`](https://github.com/janosh/matbench-discovery/blob/-/models) for example submissions. More detailed step-by-step instructions below.

### Step 1: Clone the repo

```sh
git clone https://github.com/janosh/matbench-discovery --depth 1
cd matbench-discovery
git checkout -b model-name-you-want-to-add
```

Tip: `--depth 1` only clones the latest commit, not the full `git history` which is faster if a repo contains large data files that changed over time.

### Step 2: Commit model metadata and task-specific code to the repo

Create the model family folder

```sh
mkdir -p models/<family>
```

and place the above-listed files there. The file structure should look like this:

```txt
matbench-discovery-root
└── models
    └── <family>
        ├── readme.md
        ├── <model_key>.yml
        ├── <model_key>/  # prediction artifacts
        ├── test_<model_key>_<task>.py  # only without a shared runner
        └── train_<family>.py  # optional
```

You can include arbitrary other supporting files like metadata and model features (below 10MB total to keep `git clone` time low) if they are needed to run the model or help others reproduce your results. For larger files, please upload to [Figshare](https://figshare.com) or similar and share the link in your PR description.

Discovery, phonons, molecular dynamics, and diatomics do not need per-model scripts when a model exposes an ASE calculator. Register the calculator once in [`matbench_discovery/calculators.py`](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/calculators.py), declare its isolated dependencies in `environment`, put the verified phonon protocol under `hyperparams.evaluation.kappa`, smoke-test the exact isolated command, and run the task-specific shared runner. Model-specific phonon batching or relaxation belongs in a typed adapter under `matbench_discovery/phonons/adapters`, not a forked executable.

```sh
# WBM discovery: relax/resume shards, then strictly merge and write YAML metrics
uv run models/run_discovery.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_discovery.py --print-cmd --model <model_key>)"
eval "$(uv run models/run_discovery.py --print-cmd --model <model_key> --merge-shards --write-yaml)"

# PhononDB thermal conductivity: resumable per-material records + strict merge
uv run models/run_kappa.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_kappa.py --print-cmd --model <model_key>)"
eval "$(uv run models/run_kappa.py --print-cmd --model <model_key> --merge-shards --write-yaml)"

# 20 ps NVT molecular dynamics
uv run models/run_md.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_md.py --print-cmd --model <model_key> --write-yaml)"

# homonuclear diatomic curves
uv run models/run_diatomics.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_diatomics.py --print-cmd --model <model_key> --write-yaml)"
```

For a full contiguous Slurm discovery or kappa array, the runner can infer shard selection from `SLURM_ARRAY_TASK_COUNT`, `SLURM_ARRAY_TASK_ID`, `SLURM_ARRAY_TASK_MIN`, and `SLURM_ARRAY_TASK_MAX`. Any rerun of a shard subset (contiguous or sparse) must pass the original `--n-shards`; zero-based task IDs are then reused as the original shard indices. Kappa shards are balanced deterministically by atom count and resume from atomic per-material records. Only a strict, complete `--merge-shards --write-yaml` step updates model metadata, and official YAML writes reject protocol overrides and noncanonical dataset content.

The MD runner's `--write-yaml` records model-level metrics under `metrics.md` and writes the per-system predictions to `<yyyy-mm-dd>-md-metrics.csv.gz`. Public runs compute the observable metrics only; energy/force RMSEs are maintainer-computed private-label diagnostics and are not required in external submissions. Upload generated artifacts to Figshare (or similar) and set their download URLs under the corresponding YAML `pred_file.url` fields — see [Step 3](#step-3-upload-results-files-to-figshare-or-similar) for the upload conventions.

The `Model` enum is generated during ingestion from each model YAML's required `model_key`; do not edit [`matbench_discovery/enums.py`](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/enums.py) in submission PRs. Models with `lifecycle: aborted` are omitted, while `Model.active()` selects models with `lifecycle: active`.

> [!WARNING]
> Do not include a `filter_bad_preds.py` script or equivalent. This was a historical flaw in the evaluation and now all OOM outlier prediction handling is handled internally and consistently removing the need for this filtering stage.

### Step 3: Upload results files to Figshare (or similar)

Upload the canonical prediction artifacts for each submitted task to [Figshare](https://figshare.com) (or any cloud storage with stable direct-download links) and record their URLs in the corresponding `pred_file.url`, `analysis_file.url`, `force_file.url`, or `run_info_file.url` keys. These URLs are what automated ingestion downloads, so they must resolve without authentication. During ingestion we archive copies to the project's Figshare articles and update those URLs (optionally adding `size`/`md5`).

### Step 4: Open a PR to the [Matbench Discovery repo][repo]

Commit your files to the repo on a branch called `<model_key>` and create a pull request (PR) to the Matbench repository.

```sh
git add models/<family>/<model_key>.yml models/<family>/readme.md
git commit -m 'add <model_key> to Matbench Discovery leaderboard'
```

And you're done! Once tests pass, a maintainer triggers automated ingestion on your PR (see below) and after merge your model appears on the leaderboard and in all site figures! 🎉

> [!IMPORTANT]
> Keep the **"Allow edits by maintainers"** checkbox on your PR enabled (it's on by default). The ingestion bot pushes evaluation results and updated site assets as a commit onto your PR branch — without that permission, ingestion has to be run manually by a maintainer.

### Step 5: Validate your submission

Run the ingestion script (via [`uv`](https://docs.astral.sh/uv), no extra install needed) to validate your submission before opening a PR. It runs all evaluation scripts, generates required figures, and checks the PR checklist requirements (see `scripts/ingest_model.py --help` for all flags).

#### Running the validation

```sh
uv run --with-editable . scripts/ingest_model.py <model_key>
# e.g. mace-mpa-0 or mace_mpa_0 (dashes and underscores both work)
```

This command will:

1. **Check PR requirements**: Verify model metadata, prediction URLs, shared-runner registration, and any task-specific scripts
2. **Run evaluation scripts**: Execute discovery, kappa (phonon), and diatomic metrics evaluation for your model
3. **Generate figures**: Create the required energy parity and per-element error plots

### Step 6 (optional): Copy WandB runs into our project

[Weights and Biases](https://wandb.ai) ([GitHub](https://github.com/wandb/wandb)) logs training and test runs of ML models, auto-collecting metadata like what hardware the model ran on and for how long, CPU/GPU/network utilization over that period, the exact code in the launching script, and the dependency versions installed in the environment. This helps others reproduce your results or compare computational cost, so we strongly recommend tracking all runs that went into a model submission with WandB so they can be copied over to our project at <https://wandb.ai/janosh/matbench-discovery> for everyone to inspect.

## 🤖 &thinsp; Automated ingestion (what happens after your PR is opened)

Once your submission PR looks ready, a maintainer applies the `ingest-model` label. CI ([`update-site-figs.yml`](https://github.com/janosh/matbench-discovery/blob/-/.github/workflows/update-site-figs.yml)) then runs the full pipeline — no manual work needed from you:

1. **Safe checkout** — CI checks out trusted `main` code, applies *only canonical `models/<family>/<model_key>.yml` additions, updates, renames, and removals*, and generates the `Model` enum with trusted code. Submission code is never executed in CI, which makes it safe to run this workflow on unreviewed PRs with credentials present.
2. **Evals + checklist** — trusted code validates submission metadata, downloads predictions from the YAML URLs, runs all evals (discovery, kappa, geo-opt, diatomics), and writes metrics without importing PR calculator code.
3. **Figshare archival** — only after every model passes validation, prediction + analysis files are re-uploaded to the project's own Figshare articles (one per prediction task) for longevity, and your YAML's `*_url` keys are rewritten to the archived copies.
4. **Per-model site assets** — energy/kappa parity assets (published to the GitHub release the site build downloads from), parity manifests and per-element error data are generated.
5. **Multi-model figures** — all site figure payloads (`site/src/figs/*.jsonl`) are regenerated so every page (`/tasks/discovery/tmi`, `/tasks/geo-opt`, data pages) includes your model, validated by payload shape tests before committing.
6. **One commit onto your PR** — the updated YAML, payloads and manifests are pushed to your PR branch (this is why "Allow edits by maintainers" must stay enabled), retriggering CI on the result. Merging the PR lands your model fully integrated.

Calculator changes may ship with their model YAML. Changes to payload-generating or other privileged ingestion code still require a reviewed prerequisite PR or local maintainer ingestion. For ordinary non-ingestion PRs that change the payload *format*, regenerate payloads locally and commit them in the same PR:

- Roster drift fails `tests/site/test_fig_payloads.py`'s coverage tests. After adding or updating an active model, run `uv run --with-editable . scripts/ingest_model.py <model_key> --payloads-only` to splice its freshly computed entries into the committed payloads. After removing, aborting, deprecating, or superseding a model, omit the model key and run `uv run --with-editable . scripts/ingest_model.py --payloads-only` to rebuild payloads from the full active roster.
- The `model-pr-guard` check fails any PR that changes payload-generating code without also updating committed payloads (`site/src/figs` data files or the route-local `per-element-each-errors.jsonl`): regenerate locally with `scripts/ingest_model.py --payloads-only` (no model name needed for pure format changes), or add the `payloads-unchanged` label for output-neutral refactors.
- Model-independent data-page payloads and route-local element counts are owned by `scripts/export_data_fig_payloads.py`. Run it after changing that script or `matbench_discovery/data_figs.py`; it requires the local `data/mp/2022-09-16-mp-trj-summary.json.bz2` cache. If absent, regenerate that cache from the extXYZ source with `python data/mp/eda_mp_trj.py`.

Maintainer notes: ingestion requires the repo secrets `SITE_FIGS_PAT` (classic PAT with `public_repo` scope, used to push to fork branches) and `FIGSHARE_TOKEN` (archival uploads). Re-trigger by re-applying the label, or run `uv run --with-editable . scripts/ingest_model.py <model_key> --archive` locally.

## 😵‍💫 &thinsp; Troubleshooting

Having problems? [Open an issue on GitHub](https://github.com/janosh/matbench-discovery/issues). We're happy to help! 😊

[repo]: https://github.com/janosh/matbench-discovery
