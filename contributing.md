# How to submit new models to Matbench Discovery

## 🔨 &thinsp; Installation

In an activated Python 3.14+ virtual environment, clone [the repo](https://github.com/janosh/matbench-discovery) and install `matbench-discovery`:

```zsh
git clone https://github.com/janosh/matbench-discovery --depth 1
uv pip install -e ./matbench-discovery
```

There's also a [PyPI package](https://pypi.org/project/matbench-discovery) for faster installation if you don't need the latest code changes (unlikely if you're planning to submit a model since the benchmark is under active development).

## ✨ &thinsp; How to submit a new model

To submit a new model to this benchmark and add it to our leaderboard, please create a model-submission pull request to the [`main` branch][repo] as described in [Step 4](#step-4-open-a-pr-to-the-matbench-discovery-repo). The PR must include the 3 required items:

1. You should share your model's predictions through a cloud storage service (we recommend [Figshare](https://figshare.com)) and include the stable direct-download links in your model YAML. Each artifact lives directly under `models/<family>/<model_key>/` and follows the canonical filename grammar:

   - geometry optimization: `<yyyy-mm-dd>-geo-opt.jsonl.gz` — final relaxed structures with `energy`, `converged`, `n_steps`, and WBM material IDs in [JSON Lines format](https://jsonlines.org). Optional analysis files use `<yyyy-mm-dd>-geo-opt-symprec=<symprec>-moyo=<version>.csv.gz`.
   - discovery: `<yyyy-mm-dd>-discovery.csv.gz` — compressed CSV with `material_id` and final `e_form_per_atom` columns.
   - phonons: `<yyyy-mm-dd>-phonons-kappa-103.json.gz` — predictions for exactly the 103 PhononDB material IDs. Use canonical fields `material_id`, `init_spg_num`, `final_spg_num`, `has_imag_ph_modes`, and `ph_freqs`; legacy spellings such as `mp_id`, `spg_num`, and `frequencies` are rejected rather than converted while loading. The shared runner also writes `-phonons.json.gz` (referenced as `phonon_file`): unrounded compact FC2 (primitive-atom rows, which phonopy expands losslessly) with cell bases and masses, weighted mesh spectra, seekpath bands with eigenvectors and group velocities, tetrahedron DOS, and thermal properties. Schema v3 records can be restored with `matbench_discovery.phonons.harmonic.phonopy_from_harmonic_data(record)` to evaluate new metrics without MLIP force calls. Each harmonic analysis reports its own `errors`; these also appear as `phonon_errors` in predictions and never invalidate successful conductivity. Resume skips recorded failures unless `--retry-failures` is passed. The site displays modes for successful band paths. Optional force sets and run provenance use `-forces.json.gz` and `-run-info.json`.
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

   The same shape is used for `analysis_file`, `force_file`, `phonon_file`, and `run_info_file`.

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
           (0, 0, 0),
           (0.5, 0.5, 0.5),
           (0.5, 0, 0.5),
           (0, 0.5, 0.5),
           (0.5, 0.5, 0),
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
   opt.attach(callback)  # register callback
   opt.run(fmax=0.01, steps=500)  # optimize geometry

   df_traj = pd.DataFrame(trajectory)
   df_traj.index.name = "step"
   df_traj.to_csv("trajectory.csv.gz")  # Save final structure and trajectory data
   ```

2. Reproducible inference code. If your model exposes an ASE calculator, register its constructor in [`matbench_discovery/calculators.py`](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/calculators.py) and declare its isolated `uv` dependencies in YAML `environment`. The shared runners run on Python 3.14 with `matbench-discovery`'s own `ase`, `numpy`, `pymatgen`, `phonopy`, etc., so they install only your model-specific packages from `environment.dependencies` (e.g. your package, `torch`, PyG wheels) and treat `python_version` and pins of packages `matbench-discovery` already depends on as submission provenance. Your calculator stack must therefore resolve on Python 3.14 alongside those versions; check with `--print-cmd ... --dry-run`. The shared runners [`models/run_discovery.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_discovery.py), [`models/run_kappa.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_kappa.py), [`models/run_md.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_md.py), and [`models/run_diatomics.py`](https://github.com/janosh/matbench-discovery/blob/main/models/run_diatomics.py) are the only entry points for calculator-backed models; do not add per-model forks or join scripts.

   Calculator integration belongs in the same model PR, where CI tests it.

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
   # submit as pending; ingestion sets active (later: superseded, deprecated or aborted)
   lifecycle: pending
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

   environment: # submission provenance; runners install only the model-specific deps
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
   # test_set_informed: [discovery]  # optional, see below
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

   Set `test_set_informed` to every task (`discovery`, `geo_opt`, `phonons`, `diatomics`, `md`) whose test materials or results influenced the submitted model, directly or through earlier checkpoints. That covers test materials in the training or validation data, architecture or hyperparameter choices made from test-set ablations (even on a subsample), design changes prompted by test-set errors, checkpoint selection by benchmark score, and fine-tuning judged on the test task. Explain each case in one sentence under `notes.disclosure`, citing its source. Standard training sets (MPtrj, OMat24, sAlex, ...) and choices validated only on non-test data don't count. Flagged models are marked † on the leaderboards, not excluded.

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

# 20 ps NVT molecular dynamics. Writes observable metrics only: energy/force RMSEs need
# the private labeled DynaMat references, so maintainers add them after submission
uv run models/run_md.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_md.py --print-cmd --model <model_key> --write-yaml)"

# homonuclear diatomic curves
uv run models/run_diatomics.py --print-cmd --model <model_key> --dry-run
eval "$(uv run models/run_diatomics.py --print-cmd --model <model_key> --write-yaml)"
```

For a full contiguous Slurm discovery or kappa array, the runner can infer shard selection from `SLURM_ARRAY_TASK_COUNT`, `SLURM_ARRAY_TASK_ID`, `SLURM_ARRAY_TASK_MIN`, and `SLURM_ARRAY_TASK_MAX`. Any rerun of a shard subset (contiguous or sparse) must pass the original `--n-shards`; zero-based task IDs are then reused as the original shard indices. Kappa shards are balanced deterministically by atom count and resume from atomic per-material records. Only a strict, complete `--merge-shards --write-yaml` step updates model metadata, and official YAML writes reject protocol overrides and noncanonical dataset content.

The MD runner's `--write-yaml` records model-level metrics under `metrics.md` and writes the per-system predictions to `<yyyy-mm-dd>-md-metrics.csv.gz`. Public runs compute the observable metrics only; energy/force RMSEs are maintainer-computed private-label diagnostics and are not required in external submissions. Upload generated artifacts to Figshare (or similar) and set their download URLs under the corresponding YAML `pred_file.url` fields — see [Step 3](#step-3-upload-results-files-to-figshare-or-similar) for the upload conventions.

For parallel MD runs, aggregate one selected run directory with `uv run scripts/evals/md.py --models <model_key> --md-run-dir <run_directory>`. Missing, duplicate, or unexpected systems abort the update. Without `--md-run-dir`, the evaluator only recomputes the YAML's declared prediction and preserves its download URL and checksums.

#### Inference cost provenance

Every task section (`metrics.discovery`, `metrics.diatomics`, `metrics.md`, `metrics.phonons.kappa_103` and its `harmonic` rerun) may hold one atomic `cost` mapping for the run that produced the artifact next to it. The shared runners write it on `--write-yaml`; don't hand-edit runner blocks, since recomputes re-derive them from the artifact and reject YAML costs the artifact does not back.

```yaml
cost: {source: runner, hardware: NVIDIA H200, device: cuda, dtype: float64, acceleration: [cuequivariance], n_shards: 8, n_processes: 9, run_time_sec: 540.47, load_time_sec: 96.0, max_rss_gb: 4.632, max_gpu_mem_gb: 2.729}
```

`source` is `runner` (recorded in the registered artifact: MD CSV columns, diatomics `run_metadata`, discovery/kappa `run_info_file`), `salvaged` (our earlier cluster runs predating this provenance) or `author_reported` (submitters' own hardware; if the artifact records the cost, the YAML must match it). `hardware` uses the `torch.cuda.get_device_name` spelling or `CPU (<arch>)`, `device`/`dtype` are where the weights actually ran, `acceleration` lists runner-enabled paths. `run_time_sec` is steady-state compute summed over units, `load_time_sec` model load/compile summed over the `n_processes` OS processes, memory peaks are maxed (`max_rss_gb` includes everything loaded before the timed work, for discovery mostly the WBM dataset). A new artifact replaces an `author_reported` cost only with `--replace-author-reported`. A merge over shards/systems with mixed or partial run labels (hardware, device, dtype, acceleration) drops the whole cost block with a warning.

#### Reproducibility checks

Maintainers record reruns of a registered task under `metrics.{diatomics,md}.reproducibility` with `scripts/evals/reproducibility.py --task {diatomics,md} --models <model> --run LABEL PATH [URL]`, where PATH is a canonical artifact in the model dir or a subfolder like `repeat/` (for MD a combined per-system CSV). Runs are scored with the task's own metric code against the registered `pred_file`: `rel_changes` lists each nonzero max relative change, `repeat_rel_change` (`--repeat A B`) the noise floor between identical runs. The `status` is `reproduced` when every diatomic metric, or for MD the frame-wise `energy_rmse`/`force_rmse` (chaotic trajectory metrics are only listed), changed by at most `tolerance` (fixed at 1e-3), else `drifted` (needs `--cause` and `--note`); a verdict metric missing from any file makes it `inconclusive` (needs `--note`). Registering a new `pred_file` drops the block.

The `Model` enum is built at import from every model YAML's `model_key`, so a new YAML is all a submission needs. Submit with `lifecycle: pending`: CI validates the model without requiring site payloads, and ingestion sets it `active` when it adds them. Models with `lifecycle: aborted` are omitted, while `Model.active()` selects models with `lifecycle: active`.

> [!WARNING]
> Do not include a `filter_bad_preds.py` script or equivalent. This was a historical flaw in the evaluation and now all OOM outlier prediction handling is handled internally and consistently removing the need for this filtering stage.

### Step 3: Upload results files to Figshare (or similar)

Upload the canonical prediction artifacts for each submitted task to [Figshare](https://figshare.com) (or any cloud storage with stable direct-download links) and record their URLs in the corresponding `pred_file.url`, `analysis_file.url`, `force_file.url`, `phonon_file.url`, or `run_info_file.url` keys. These URLs are what automated ingestion downloads, so they must resolve without authentication. During ingestion we archive copies to the project's Figshare articles and update those URLs (optionally adding `size`/`md5`).

### Step 4: Open a PR to the [Matbench Discovery repo][repo]

Commit and push your files on a branch called `<model_key>`:

```sh
git add models/<family>/<model_key>.yml models/<family>/readme.md
git commit -m 'add <model_key> to Matbench Discovery leaderboard'
git push -u origin <model_key>
```

Then replace `YOUR_GITHUB_USER` and `MODEL_KEY` in this branch-specific URL. The
explicit base and head branches ensure GitHub loads the model checklist:

```text
https://github.com/janosh/matbench-discovery/compare/main...YOUR_GITHUB_USER:MODEL_KEY?expand=1&template=new_model.md
```

Alternatively, use the GitHub CLI from the repository checkout:

```sh
gh pr create --repo janosh/matbench-discovery --base main \
  --head YOUR_GITHUB_USER:MODEL_KEY --draft \
  --template .github/PULL_REQUEST_TEMPLATE/new_model.md
```

And you're done! Once tests pass, a maintainer reviews and ingests your model (see below), and after merge it appears on the leaderboard and in all site figures! 🎉

> [!IMPORTANT]
> Keep the **"Allow edits by maintainers"** checkbox on your PR enabled (it's on by default). The maintainer pushes evaluation results and updated site assets as a commit onto your PR branch.

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

## 🤖 &thinsp; Ingestion (what happens after your PR is opened)

After review, a maintainer runs `uv run --with-editable . scripts/ingest_model.py <model_key> --archive` and pushes the result onto your PR as one commit. The script validates the submission, runs all evals, archives your prediction files to the project's Figshare articles (rewriting your YAML's URLs to the archived copies), publishes the energy/kappa parity assets to the GitHub release, regenerates every site figure payload and sets `lifecycle: active`. Merging the PR lands your model fully integrated.

For PRs that change the payload *format* or the model roster, regenerate payloads locally and commit them in the same PR:

- Roster drift fails `tests/site/test_fig_payloads.py`'s coverage tests. After adding or updating an active model, run `uv run --with-editable . scripts/ingest_model.py <model_key> --payloads-only` to splice its freshly computed entries into the committed payloads. After removing, aborting, deprecating, or superseding a model, omit the model key and run `uv run --with-editable . scripts/ingest_model.py --payloads-only --full-roster` to rebuild payloads from the full active roster.
- Payloads record hashes of the Python sources that computed them, ignoring comments, docstrings and formatting. A code change to one of those files fails `tests/site/test_fig_payloads.py` until you regenerate with `uv run --with-editable . scripts/ingest_model.py --payloads-only --full-roster`.
- Parity generators are separate because their large assets live on the GitHub release. Run `uv run --with-editable . site/scripts/generate-energy-parity-assets.py` and `uv run --with-editable ".[phonons]" site/scripts/generate-kappa-parity-assets.py`, publish the generated files with `uv run --with-editable . scripts/ingest_model.py --publish-parity`, then commit both manifests.
- Model-independent data-page payloads and route-local element counts are owned by `scripts/export_data_fig_payloads.py`. Run it after changing that script or `matbench_discovery/data_figs.py`; it requires the local `data/mp/2022-09-16-mp-trj-summary.json.bz2` cache. If absent, regenerate that cache from the extXYZ source with `uv run data/mp/eda_mp_trj.py`.

Maintainer note: `--archive` needs a Figshare API token in `$FIGSHARE_TOKEN` or as a `figshare_token=` line in `site/.env`, and `--publish-parity` needs `gh` with write access to this repo's releases.

## 😵‍💫 &thinsp; Troubleshooting

Having problems? [Open an issue on GitHub](https://github.com/janosh/matbench-discovery/issues). We're happy to help! 😊

[repo]: https://github.com/janosh/matbench-discovery
