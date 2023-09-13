<script>
  import { name, repository as repo, homepage } from "$site/package.json"
  import figshare_urls from "$root/data/figshare/1.0.0.json"
  import { Tooltip } from 'svelte-zoo'

  const ppd_doc_url = `https://github.com/materialsproject/pymatgen/blob/v2023.5.10/pymatgen/analysis/phase_diagram.py#L1480-L1814`
  const ppd_link = `<a href=${ppd_doc_url}>PatchedPhaseDiagram</a>`
  const cse_doc_url = `https://github.com/materialsproject/pymatgen/blob/v2023.5.10/pymatgen/entries/computed_entries.py#L579-L722`
  const cse_link = `<a href=${cse_doc_url}>ComputedStructureEntry</a>`

  const descriptions = {
    alignn_checkpoint: "ALIGNN model trained on <code>mp_computed_structure_entries</code>",
    mace_checkpoint: "MACE model trained on the MPtrj dataset (https://figshare.com/articles/dataset/23713842)",
    mp_computed_structure_entries:
      `JSON-Serialized MP ${cse_link} objects containing relaxed structures and DFT final energies`,
    mp_elemental_ref_entries: `Minimum energy ComputedEntry for each element in MP`,
    mp_energies: `Materials Project formation energies and energies above convex hull`,
    mp_patched_phase_diagram:
      `${ppd_link} constructed from all MP ComputedStructureEntries`,
    wbm_computed_structure_entries: `Materials Project computed structure entries`,
    wbm_initial_structures: `Unrelaxed WBM structures`,
    wbm_cses_plus_init_structs: `Both unrelaxed and DFT-relaxed WBM structures, the latter stored with their final VASP energies as ${cse_link}`,
    wbm_summary:
      `Computed material properties only, no structures. Available properties are VASP energy, formation energy, energy above the convex hull, volume, band gap, number of sites per unit cell, and more.`,
  }
  const desc_keys = Object.keys(descriptions).sort()
  const figshare_keys = Object.keys(figshare_urls).sort()
  const missing = figshare_keys.filter((key) => !desc_keys.includes(key))
  if (missing.length > 0) {
    throw `descriptions must contain all figshare_urls keys, missing=${missing}`
  }
</script>

# How to contribute

## 🔨 &thinsp; Installation

To download the training and test sets for this benchmark, we recommend installing our [PyPI package](https://pypi.org/project/{name}):

```zsh
pip install matbench-discovery
```

## 📙 &thinsp; Usage

When you access an attribute of the `DATA_FILES` class, it automatically downloads and caches the corresponding data file. For example:

```py
from matbench_discovery.data import DATA_FILES, load

# available data files
assert sorted(DATA_FILES) == [
    "mp_computed_structure_entries",
    "mp_elemental_ref_entries",
    "mp_energies",
    "mp_patched_phase_diagram",
    "wbm_computed_structure_entries",
    "wbm_cses_plus_init_structs",
    "wbm_initial_structures",
    "wbm_summary",
]

# version defaults to latest, set a specific version to avoid breaking changes
df_wbm = load("wbm_summary", version="1.0.0")

assert df_wbm.shape == (256963, 15)

assert list(df_wbm) == [
    "formula",
    "n_sites",
    "volume",
    "uncorrected_energy",
    "e_form_per_atom_wbm",
    "e_above_hull_wbm",
    "bandgap_pbe",
    "wyckoff_spglib",
    "uncorrected_energy_from_cse",
    "e_correction_per_atom_mp2020",
    "e_correction_per_atom_mp_legacy",
    "e_form_per_atom_uncorrected",
    "e_form_per_atom_mp2020_corrected",
    "e_above_hull_mp2020_corrected_ppd_mp",
    "site_stats_fingerprint_init_final_norm_diff",
]
```

`"wbm-summary"` columns:

1. **`formula`**: A compound's unreduced alphabetical formula
1. **`n_sites`**: Number of sites in the structure's unit cell
1. **`volume`**: Relaxed structure volume in cubic Angstrom
1. **`uncorrected_energy`**: Raw VASP-computed energy
1. **`e_form_per_atom_wbm`**: Original formation energy per atom from [WBM paper]
1. **`e_above_hull_wbm`**: Original energy above the convex hull in (eV/atom) from [WBM paper]
1. **`wyckoff_spglib`**: Aflow label strings built from spacegroup and Wyckoff positions of the relaxed structure as computed by [spglib](https://spglib.readthedocs.io/en/latest/python-spglib.html?highlight=get_symmetry_dataset#get-symmetry-dataset).
1. **`bandgap_pbe`**: PBE-level DFT band gap from [WBM paper]
1. **`uncorrected_energy_from_cse`**: Uncorrected DFT energy stored in `ComputedStructureEntries`. Should be the same as `uncorrected_energy`. There are 2 cases where the absolute difference reported in the summary file and in the computed structure entries exceeds 0.1 eV (`wbm-2-3218`, `wbm-1-56320`) which we attribute to rounding errors.
1. **`e_form_per_atom_uncorrected`**: Uncorrected DFT formation energy per atom in eV/atom.
1. **`e_form_per_atom_mp2020_corrected`**: Matbench Discovery takes these as ground truth for the formation energy. The result of applying the [MP2020 energy corrections](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility) (latest correction scheme at time of release) to `e_form_per_atom_uncorrected`.
1. **`e_correction_per_atom_mp2020`**: [`MaterialsProject2020Compatibility`](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility) energy corrections in eV/atom.
1. **`e_correction_per_atom_mp_legacy`**: Legacy [`MaterialsProjectCompatibility`](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProjectCompatibility) energy corrections in eV/atom. Having both old and new corrections allows updating predictions from older models like MEGNet that were trained on MP formation energies treated with the old correction scheme.
1. **`e_above_hull_mp2020_corrected_ppd_mp`**: Energy above hull distances in eV/atom after applying the MP2020 correction scheme. The convex hull in question is the one spanned by all ~145k Materials Project `ComputedStructureEntries`. Matbench Discovery takes these as ground truth for material stability. Any value above 0 is assumed to be an unstable/metastable material.
1. **`site_stats_fingerprint_init_final_norm_diff`**: The norm of the difference between the initial and final site fingerprints. This is a volume-independent measure of how much the structure changed during DFT relaxation. Uses the `matminer` [`SiteStatsFingerprint`](https://github.com/hackingmaterials/matminer/blob/33bf112009b67b108f1008b8cc7398061b3e6db2/matminer/featurizers/structure/sites.py#L21-L33) (v0.8.0).

## 📥 &thinsp; Direct Download

You can download the data files from Figshare:

<ol>
  {#each Object.entries(figshare_urls) as [key, lst]}
    {@const [href, file_name] = lst}
    <li>
      <Tooltip text={file_name}>
        <a {href}>{key}</a>:
      </Tooltip>
      {@html descriptions[key]}
    </li>
  {/each}
</ol>

[wbm paper]: https://nature.com/articles/s41524-020-00481-6

## ✨ &thinsp; How to submit a new model

To submit a new model to this benchmark and add it to our leaderboard, please create a pull request to the [`main` branch]({repo}) that includes at least these 3 required files:

1. `<yyyy-mm-dd>-<model_name>-preds.(json|csv).gz`: Your model's energy predictions for all ~250k WBM compounds as compressed JSON or CSV. The recommended way to create this file is with `pandas.DataFrame.to_{json|csv}('<yyyy-mm-dd>-<model_name>-preds.(json|csv).gz')`. JSON is preferred over CSV if your model not only predicts energies (floats) but also objects like relaxed structures. See e.g. [M3GNet]({repo}/blob/-/models/m3gnet/test_m3gnet.py) and [CHGNet]({repo}/blob/-/models/chgnet/test_chgnet.py) test scripts.
1. `test_<model_name>.(py|ipynb)`: The Python script or Jupyter notebook that generated the energy predictions. Ideally, this file should have comments explaining at a high level what the code is doing and how the model works so others can understand and reproduce your results. If the model deployed on this benchmark was trained specifically for this purpose (i.e. if you wrote any training/fine-tuning code while preparing your PR), please also include it as `train_<model_name>.(py|ipynb)`.
1. `metadata.yml`: A file to record all relevant metadata of your algorithm like model name and version, authors, package requirements, links to publications, notes, etc. Here's a template:

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

   training_set:
     title: MPtrj
     url: https://figshare.com/articles/dataset/23713842
     size: 1_580_395

   notes: # notes can have any key, be multiline and support markdown.
     description: This is how my model works...
     steps: |
      Optional *free-form* [markdown](example.com) notes.
   ```

   Arbitrary other keys can be added as needed.

Please see any of the subdirectories in [`models/`]({repo}/blob/-/models) for example submissions. More detailed step-by-step instructions below.

### Step 1: Clone the repo

```sh
git clone https://github.com/janosh/matbench-discovery --depth 1
cd matbench-discovery
git checkout -b model-name-you-want-to-add
```

Use `--depth 1` to only clone the latest commit and save bandwidth. No need to load previous versions of large data files.

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

This information can be useful for others looking to reproduce your results or compare their model to yours i.t.o. computational cost. We therefore strongly recommend tracking all runs that went into a model submission with WandB so that the runs can be copied over to our WandB project at <https://wandb.ai/janosh/matbench-discovery> for everyone to inspect. This also allows us to include your model in more detailed analysis (see [SI]({homepage}/preprint#supplementary-information).

## 😵‍💫 &thinsp; Troubleshooting

Having problems? Please [open an issue on GitHub]({repo}/issues). We're happy to help! 😊
