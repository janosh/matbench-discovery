<script>
  import { data_files } from '$lib'
</script>

## 📙 &thinsp; Downloading Data Files

The `DataFiles` class in `matbench_discovery.data` provides automatically downloads and locally caches data files needed for different model training and testing tasks. For example the following code snippet shows how to download and access the WBM summary data:

```py
from matbench_discovery.data import DataFiles, ase_atoms_from_zip
import pandas as pd

df_wbm = pd.read_csv(DataFiles.wbm_summary.path)
# confirm test set size
assert df_wbm.shape == (256_963, 18)
# available columns in WBM summary data
assert tuple(df_wbm) == (
    "material_id",
    "formula",
    "n_sites",
    "volume",
    "uncorrected_energy",
    "e_form_per_atom_wbm",
    "e_above_hull_wbm",
    "bandgap_pbe",
    "wyckoff_spglib_initial_structure",
    "uncorrected_energy_from_cse",
    "e_correction_per_atom_mp2020",
    "e_correction_per_atom_mp_legacy",
    "e_form_per_atom_uncorrected",
    "e_form_per_atom_mp2020_corrected",
    "e_above_hull_mp2020_corrected_ppd_mp",
    "site_stats_fingerprint_init_final_norm_diff",
    "wyckoff_spglib",
    "unique_prototype"
)

# WBM initial structures in pymatgen JSON format
df_init_structs = pd.read_json(DataFiles.wbm_initial_structures.path, lines=True)
assert tuple(df_init_structs) == ("material_id", "formula_from_cse", "initial_structure")
# WBM initial structures as ASE Atoms
wbm_init_atoms = ase_atoms_from_zip(DataFiles.wbm_initial_atoms.path)
assert len(wbm_init_atoms) == 256_963
```

`"wbm-summary"` columns:

1. **`formula`**: A compound's unreduced alphabetical formula
1. **`n_sites`**: Number of sites in the structure's unit cell
1. **`volume`**: Relaxed structure volume in cubic Angstrom
1. **`uncorrected_energy`**: Raw VASP-computed energy
1. **`e_form_per_atom_wbm`**: Original formation energy per atom from [WBM paper]
1. **`e_above_hull_wbm`**: Original energy above the convex hull in (eV/atom) from [WBM paper]
1. **`wyckoff_spglib`**: Aflow label strings built from spacegroup and Wyckoff positions of the DFT-relaxed structure as computed by [spglib](https://spglib.readthedocs.io/en/stable/api/python-api.html#spglib.spglib.get_symmetry_dataset).
1. **`wyckoff_spglib_initial_structure`**: Same as `wyckoff_spglib` but computed from the initial structure.
1. **`bandgap_pbe`**: PBE-level DFT band gap from [WBM paper]
1. **`uncorrected_energy_from_cse`**: Uncorrected DFT energy stored in `ComputedStructureEntries`. Should be the same as `uncorrected_energy`. There are 2 cases where the absolute difference reported in the summary file and in the computed structure entries exceeds 0.1 eV (`wbm-2-3218`, `wbm-1-56320`) which we attribute to rounding errors.
1. **`e_form_per_atom_uncorrected`**: Uncorrected DFT formation energy per atom in eV/atom.
1. **`e_form_per_atom_mp2020_corrected`**: Matbench Discovery takes these as ground truth for the formation energy. The result of applying the [MP2020 energy corrections][`MaterialsProject2020Compatibility`] (latest correction scheme at time of release) to `e_form_per_atom_uncorrected`.
1. **`e_correction_per_atom_mp2020`**: [`MaterialsProject2020Compatibility`] energy corrections in eV/atom.
1. **`e_correction_per_atom_mp_legacy`**: Legacy [`MaterialsProjectCompatibility`] energy corrections in eV/atom. Having both old and new corrections allows updating predictions from older models like MEGNet that were trained on MP formation energies treated with the old correction scheme.
1. **`e_above_hull_mp2020_corrected_ppd_mp`**: Energy above hull distances in eV/atom after applying the MP2020 correction scheme. The convex hull in question is the one spanned by all ~145k Materials Project `ComputedStructureEntries`. Matbench Discovery takes these as ground truth for material stability. Any value above 0 is assumed to be an unstable/metastable material.
1. **`site_stats_fingerprint_init_final_norm_diff`**: The norm of the difference between the initial and final site fingerprints. This is a volume-independent measure of how much the structure changed during DFT relaxation. Uses the `matminer` [`SiteStatsFingerprint`](https://github.com/hackingmaterials/matminer/blob/33bf1120/matminer/featurizers/structure/sites.py#L21-L33) (v0.8.0).

[`MaterialsProject2020Compatibility`]: https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/compatibility.py#L823
[`MaterialsProjectCompatibility`]: https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/compatibility.py#L766

## 📥 &thinsp; Direct Download

You can also directly download Matbench Discovery data files from [this Figshare article](https://figshare.com/articles/dataset/23713842).

<ol class="data-files-list">
{#each Object.entries(data_files).filter(([key]) => !key.startsWith(`_`)) as [key, { url, path, html }]}
    <li style="margin-top: 1ex;">
    <strong><code>{key}</code></strong>
    (<a href={url}>{path}</a>)<br />
    {@html html}
    </li>
{/each}
</ol>

To train an interatomic potential, we recommend the [**MPtrj dataset**](https://figshare.com/articles/dataset/23713842) which was created to train [CHGNet](https://www.nature.com/articles/s42256-023-00716-3). With thanks to [Bowen Deng](https://scholar.google.com/citations?user=PRPXA0QAAAAJ) for cleaning and releasing this dataset. It was created from the [2021.11.10](https://docs.materialsproject.org/changes/database-versions#v2021.11.10) release of Materials Project and therefore constitutes a slightly smaller but valid subset of the allowed [2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28) MP release that is our training set.

[wbm paper]: https://nature.com/articles/s41524-020-00481-6
