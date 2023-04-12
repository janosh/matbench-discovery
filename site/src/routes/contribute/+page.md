<script>
  import { name, repository as repo, homepage } from "$site/package.json";
</script>

# How to contribute

## 🔨 &thinsp; Installation

The recommended way to acquire the training and test sets for this benchmark is through its Python package [available on PyPI](https://pypi.org/project/{name}):

```zsh
pip install matbench-discovery
```

## 📙 &thinsp; Usage

This example script downloads the training and test data for training a model:

```py
from matbench_discovery.data import load_train_test
from matbench_discovery.data import df_wbm, DATA_FILES

# any subset of these keys can be passed to load_train_test()
assert sorted(DATA_FILES) == [
    "mp-computed-structure-entries",
    "mp-elemental-ref-energies",
    "mp-energies",
    "mp-patched-phase-diagram",
    "wbm-computed-structure-entries",
    "wbm-initial-structures",
    "wbm-summary",
]

df_wbm = load_train_test("wbm-summary", version="v1.0.0")

assert df_wbm.shape == (256963, 15)

assert list(df_wbm) == [
    "formula",
    "n_sites",
    "volume",
    "uncorrected_energy",
    "e_form_per_atom_wbm",
    "e_above_hull_wbm",
    "bandgap_pbe",
    "uncorrected_energy_from_cse",
    "e_correction_per_atom_mp2020",
    "e_correction_per_atom_mp_legacy",
    "e_above_hull_mp2020_corrected_ppd_mp",
    "e_form_per_atom_uncorrected",
    "e_form_per_atom_mp2020_corrected",
    "wyckoff_spglib",
]
```

`"wbm-summary"` column glossary:

1. `formula`: A compound's unreduced alphabetical formula
1. `n_sites`: Number of sites in the structure's unit cell
1. `volume`: Relaxed structure volume in cubic Angstrom
1. `uncorrected_energy`: Raw VASP-computed energy
1. `e_form_per_atom_wbm`: Original formation energy per atom from [WBM paper]
1. `e_hull_wbm`: Original energy above the convex hull in (eV/atom) from [WBM paper]
1. `bandgap_pbe`: PBE-level DFT band gap from [WBM paper]
1. `uncorrected_energy_from_cse`: Should be the same as `uncorrected_energy`. There are 2 cases where the absolute difference reported in the summary file and in the computed structure entries exceeds 0.1 eV (`wbm-2-3218`, `wbm-1-56320`) which we attribute to rounding errors.
1. `e_form_per_atom_mp2020_corrected`: Matbench Discovery takes these as ground truth for the formation energy. Includes MP2020 energy corrections (latest correction scheme at time of release).
1. `e_correction_per_atom_mp2020`: [`MaterialsProject2020Compatibility`](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility) energy corrections in eV/atom.
1. `e_correction_per_atom_mp_legacy`: Legacy [`MaterialsProjectCompatibility`](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProjectCompatibility) energy corrections in eV/atom. Having both old and new corrections allows updating predictions from older models like MEGNet that were trained on MP formation energies treated with the old correction scheme.
1. `e_above_hull_mp2020_corrected_ppd_mp`: Energy above hull distances in eV/atom after applying the MP2020 correction scheme. The convex hull in question is the one spanned by all ~145k Materials Project `ComputedStructureEntries`. Matbench Discovery takes these as ground truth for material stability. Any value above 0 is assumed to be an unstable/metastable material.
<!-- TODO document remaining columns, or maybe drop them from df -->

## 📥 &thinsp; Direct Download

You can also download the data files directly from GitHub:

1. [`2022-10-19-wbm-summary.csv`]({repo}/blob/-/data/wbm/2022-10-19-wbm-summary.csv): Computed material properties only, no structures. Available properties are VASP energy, formation energy, energy above the convex hull, volume, band gap, number of sites per unit cell, and more. e_form_per_atom and e_above_hull each have 3 separate columns for old, new and no Materials
1. [`2022-10-19-wbm-init-structs.json`]({repo}/blob/-/data/wbm/2022-10-19-wbm-init-structs.json): Unrelaxed WBM structures
1. [`2022-10-19-wbm-cses.json`]({repo}/blob/-/data/wbm/2022-10-19-wbm-cses.json): Relaxed WBM structures along with final VASP energies
1. [`2023-01-10-mp-energies.json.gz`]({repo}/blob/-/data/mp/2023-01-10-mp-energies.json.gz): Materials Project formation energies and energies above convex hull
1. [`2023-02-07-mp-computed-structure-entries.json.gz`]({repo}/blob/-/data/mp/2023-02-07-mp-computed-structure-entries.json.gz): Materials Project computed structure entries
1. [`2023-02-07-ppd-mp.pkl.gz`]({repo}/blob/-/data/mp/2023-02-07-ppd-mp.pkl.gz): [PatchedPhaseDiagram](https://pymatgen.org/pymatgen.analysis.phase_diagram.html#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram) constructed from all MP ComputedStructureEntries
1. [`2022-09-19-mp-elemental-reference-entries.json`]({repo}/blob/-/data/mp/2022-09-19-mp-elemental-reference-entries.json): Minimum energy PDEntries for each element present in the Materials Project

[wbm paper]: https://nature.com/articles/s41524-020-00481-6

## ✨ &thinsp; How to submit a new model

To deploy a new model on this benchmark and add it to our leaderboard, please create a pull request to the `main` branch of [{repo}]({repo}) that includes at least these 3 required files:

1. `<yyyy-mm-dd>-<model_name>-preds.(json|csv).gz`: Your model's energy predictions for all ~250k WBM compounds as compressed JSON or CSV. The recommended way to create this file is with `pandas.DataFrame.to_{json|csv}('<yyyy-mm-dd>-<model_name>-preds.(json|csv).gz')`. JSON is preferred over CSV if your model not only predicts energies (floats) but also Python objects like e.g. pseudo-relaxed structures (see the M3GNet and BOWSR test scripts).
1. `test_<model_name>.(py|ipynb)`: The Python script or Jupyter notebook that generated the energy predictions. Ideally, this file should have comments explaining at a high level what the code is doing and how the model works so others can understand and reproduce your results. If the model deployed on this benchmark was trained specifically for this purpose (i.e. if you wrote any training/fine-tuning code while preparing your PR), please also include it as `train_<model_name>.(py|ipynb)`.
1. `metadata.yml`: A file to record all relevant metadata of your algorithm like model name and version, authors, package requirements, relevant citations/links to publications, notes, etc. Here's a template:

   ```yml
   # metadata.yml template
   model_name: My fancy model # required
   model_version: 1.0.0 # required
   matbench_discovery_version: 1.0 # required
   date_added: "2023-01-01" # required
   authors: # required (only name, other keys are optional)
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
   repo: https://github.com/<user>/<repo> # required
   url: https://<model-docs-or-similar>.org
   doi: https://doi.org/10.5281/zenodo.0000000
   preprint: https://arxiv.org/abs/xxxx.xxxxx
   requirements: # strongly recommended
     torch: 1.13.0
     torch-geometric: 2.0.9
     ...
   notes: # notes can have any key, be multiline and support markdown.
     description: This is how my model works...
     steps: |
      Optional free form multi-line notes that can help others reproduce your results.
   ```

   Arbitrary other keys can be added as needed.

Please see any of the subdirectories in [`models/`]({repo}/tree/main/models) for example submissions. More detailed step-by-step instructions below.

### Step 1: Clone the repo

```sh
git clone https://github.com/janosh/matbench-discovery
cd matbench-discovery
git checkout -b model-name-you-want-to-add
```

### Step 2: Commit model preds, script and metadata

Create a new folder

```sh
mkdir models/<model_name>
```

and place the above-listed files there. The file structure should look like this:

```txt
matbench-discovery-root
└── models
    └── <model_name>
        ├── metadata.yml
        ├── <yyyy-mm-dd>-<model_name>-preds.(json|csv).gz
        ├── test_<model_name>.py
        ├── readme.md  # optional
        └── train_<model_name>.py  # optional
```

You can include arbitrary other supporting files like metadata and model features (below 10MB to keep `git clone` time low) if they are needed to run the model or help others reproduce your results. For larger files, please upload to [Figshare](https://figshare.com) or similar and link them somewhere in your files.

### Step 3: Open a PR to the [Matbench Discovery repo]({repo})

Commit your files to the repo on a branch called `<model_name>` and create a pull request (PR) to the Matbench repository.

```sh
git add -a models/<model_name>
git commit -m 'add <model_name> to Matbench Discovery leaderboard'
```

And you're done! Once tests pass and the PR is merged, your model will be added to the leaderboard! 🎉

### Step 4 (optional): Copy WandB runs into our project

[Weights and Biases](https://wandb.ai) is a great tool for logging training and test runs of ML models. It's free, (partly) [open source](https://github.com/wandb/wandb) and offers a [special plan for academics](https://wandb.ai/site/research). It auto-collects metadata like

- what hardware the model is running on
- and for how long,
- what the CPU, GPU and network utilization was over that period,
- the exact code in the script that launched the run, and
- which versions of dependencies were installed in the environment your model ran in.

This information can be very useful for someone looking to reproduce your results or compare their model to yours i.t.o. computational cost. We therefore strongly recommend tracking all runs that went into a model submission to Matbench Discovery with WandB so that the runs can be copied over to our WandB project at <https://wandb.ai/janosh/matbench-discovery> for everyone to inspect. This also allows us to include your model in more detailed analysis found in the [SI]({homepage}/si).

## 😵‍💫 &thinsp; Troubleshooting

Having problems using or contributing to the project? Please [open an issue on GitHub]({repo}/issues) or [start a discussion]({repo}/discussions) for open-ended conversations. Happy to help! 😊
