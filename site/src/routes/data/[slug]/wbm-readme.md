<script lang="ts">
  import type { Snippet } from 'svelte'

  // snippets the page passes in, rendered below where each belongs
  let {
    hist_e_form_per_atom,
    hist_wbm_hull_dist,
    wbm_elements_heatmap,
    mp_elements_heatmap,
    mp_trj_elements_heatmap,
    spacegroup_sunbursts,
  }: Record<string, Snippet | undefined> = $props()
</script>

<!-- markdownlint-disable MD051 -->

## WBM Benchmark Details

[The **WBM dataset**](/data/wbm) was published in [Predicting stable crystalline compounds using chemical similarity][wbm paper] (nat comp mat, Jan 2021). The authors generated 257,487 structures through single-element substitutions on Materials Project (MP) source structures, choosing the replacement element by chemical similarity from a matrix data-mined from the [Inorganic Crystal Structure Database (ICSD)](https://icsd.products.fiz-karlsruhe.de).

The novel structures were relaxed with MP-compatible VASP inputs (`pymatgen`'s [`MPRelaxSet`](https://github.com/materialsproject/pymatgen/blob/c4998d92/pymatgen/io/vasp/MPRelaxSet.yaml)) and identical POTCARs, aiming for a database of MP-compatible novel crystals. Any degradation in model performance from training to test set should therefore largely result from extrapolation error rather than covariate shift in the underlying data.

The authors performed 5 rounds of element substitution starting from structures on the MP convex hull, each time relaxing all generated structures and adding those on the convex hull back to the pool of parent structures for the next round. The paper reports 18,479 of its 189,981 structures (9.7%) on the convex hull of MP and WBM structures combined.

Since repeated substitutions should on average increase chemical dissimilarity, these 5 iterations are a unique and compelling feature: they allow testing on successively more out-of-distribution (OOD) data slices, i.e. checking how model performance degrades on structures increasingly dissimilar from the training set (restricted to the [MP v2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28) database release or earlier for all models in this benchmark).

## 🆔 &thinsp; About the IDs

The first integer in each material ID (1 to 5, right after the prefix `wbm-`) is the substitution count, i.e. how many times one element was replaced with a chemically similar one starting from an MP source structure. The higher this number, the more chance to diffuse away from the MP training set into uncharted regions of materials space. On average, we expect model performance to decrease with substitution step count as structures further from the training set require more extrapolation.

The 2nd integer counts the materials within each iteration. It is not always consecutive since a small number of materials (~0.2%) were removed by the data-cleaning steps below, so don't be surprised to find an ID like `wbm-3-5674` followed by `wbm-3-5678`.

## 🪓 &thinsp; Data Processing Steps

The full set of processing steps used to curate the WBM test set from the raw data files ([URLs listed below](#links-to-wbm-files)) is in [`scripts/data/compile_wbm_test_set.py`](https://github.com/janosh/matbench-discovery/blob/main/scripts/data/compile_wbm_test_set.py):

- re-format material IDs: `step_1-0->wbm-1-1`, `step_1-1->wbm-1-2`, ...
- remove 6 pathological structures inserted towards the end of step 3 (`step_3_70802`, `_70803`, `_70825`, `_70826`, `_70828`, `_70829`): present in the summary and initial-structure files but not in the DFT-relaxed [`ComputedStructureEntries`](https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/computed_entries.py#L536), with 0 volume, energy, etc. in the summary. Dropping them aligns initial structures to their `ComputedStructureEntries`
- remove formation energy outliers below -5 and above 5 eV/atom (502 and 22 of 257,487 crystals respectively, including an anomaly of 500 structures at exactly -10 eV/atom)

  > WBM Formation energy distribution. 524 materials outside dashed lines were discarded.

  {@render hist_e_form_per_atom?.()}

- apply the [`MaterialsProject2020Compatibility`](https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/compatibility.py#L823) energy correction scheme to the formation energies
- compute energy to the Materials Project convex hull constructed from all MP `ComputedStructureEntries` queried on 2023-02-07 ([database release 2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28))

> Correction note: the `e_form_per_atom_mp2020_corrected` and `e_above_hull_mp2020_corrected_ppd_mp` columns track the current `pymatgen` `MaterialsProject2020Compatibility` defaults, applied to both the WBM entries and the MP entries spanning the convex hull, matching the corrections applied to model predictions. Since the 2026-11-19 summary this includes [pymatgen#3803](https://github.com/materialsproject/pymatgen/pull/3803), which only assigns an anion correction to elements with oxidation state <= -1. Relative to the 2023-12-13 summary, that changed the corrections of 3,044 WBM entries (up to 0.49 eV/atom) and 329 MP hull entries and flipped the stability label of 605 WBM entries. See [#358](https://github.com/janosh/matbench-discovery/issues/358) for the history.

Run `uv run --with-editable . scripts/data/compile_wbm_test_set.py` to download the raw WBM structures and rebuild the summary. The compiler recomputes energy corrections, formation energies (from the original four-decimal elemental energies), MP hull distances (`pymatgen` `PatchedPhaseDiagram` over the MP2020-corrected MP entries, with an exact `PhaseDiagram` for chemical systems no patch covers) and unique-prototype membership. It reuses the two spglib labels and Matminer fingerprint distance from the 2023 release ([2023 WBM reference](https://figshare.com/files/64706751)) and verifies the resulting cohort against every published material ID. If you find

- any questionable structures or data records in the released test set, or
- inconsistencies between the files on Figshare vs the output of that script,

please [raise an issue](https://github.com/janosh/matbench-discovery/issues).

Number of materials per element-substitution iteration before and after processing:

| step       | 1      | 2      | 3      | 4      | 5      | total   |
| ---------- | ------ | ------ | ------ | ------ | ------ | ------- |
| **before** | 61,848 | 52,800 | 79,205 | 40,328 | 23,308 | 257,489 |
| **after**  | 61,466 | 52,755 | 79,160 | 40,314 | 23,268 | 256,963 |

**before** counts `ComputedStructureEntries`. Two step-5 entries (`wbm-5-23166`, `wbm-5-23294`) have no initial structure, leaving the 257,487 materials the processing starts from.

### Prototype Analysis

A common issue in high-throughput screening is deciding when two materials output by a workflow should be treated as the same. The crudest option, comparing reduced formulae, is too coarse-grained to be useful in many cases: two materials can share a formula but have different crystal structures, the canonical example being graphite and diamond. Alternatively, materials can be matched by fingerprinting their crystal structures and determining their similarity through some choice of kernel function and cutoff, but this can be expensive to compute for large datasets and requires selecting multiple hyperparameters. We instead determine each crystal's prototype from the Wyckoff positions of its elements, a middle ground between the two in cost, number of hyperparameters and ability to distinguish polymorphs.

The published cohort uses AFLOW-style labels computed with Aviary and spglib ([`get_protostructure_label`], whose moyopy-based successor is `get_protostructure_label_from_moyopy` in `pymatgen.analysis.prototypes`). The compiler reads these historical labels from the 2023 release instead of recomputing them, since newer symmetry-detection versions can assign different space groups and Wyckoff orbits to distorted structures. [Original implementation][original-prototypes]. Having determined the prototypes for both MP and WBM, we first removed any WBM structure sharing a prototype with any MP structure, then dropped all but the lowest-energy structure for each remaining unique WBM prototype label. WBM prototypes were determined from the relaxed structures to avoid inflating the metrics with duplicates from different initial prototypes that fall into the same basin of attraction during relaxation. In total, this approach removes 41,475 candidates from the WBM dataset.

Number of materials per element-substitution iteration before and after prototype filtering:

| step       | 1      | 2      | 3      | 4      | 5      | total   |
| ---------- | ------ | ------ | ------ | ------ | ------ | ------- |
| **before** | 61,466 | 52,755 | 79,160 | 40,314 | 23,268 | 256,963 |
| **after**  | 54,209 | 45,979 | 66,528 | 34,531 | 14,241 | 215,488 |

[`get_protostructure_label`]: https://github.com/janosh/matbench-discovery/blob/26d3a19073d/matbench_discovery/structure/prototype.py#L104-L193
[original-prototypes]: https://github.com/janosh/matbench-discovery/blob/a7f693d1/data/wbm/compile_wbm_test_set.py#L587-L616

## 🔗 &thinsp; Links to WBM Files

Links to raw WBM data files have proliferated. This section tracks them.

Initial structures (after element substitution, before DFT relaxation) were sent as Google Drive links via email by Hai-Chen Wang on 2021-09-01.

### Google Drive

| Google Drive links | [step 1](https://drive.google.com/file/d/1ZUgtYwrfZn_P8bULWRtTXepyAxHVxS5C) | [step 2](https://drive.google.com/file/d/1-3uu2AcARJxH7GReteGVASZTuttFGiW_) | [step 3](https://drive.google.com/file/d/1hc5BvDiFfTu_tc5F8m7ONSw2OgL9vN6o) | [step 4](https://drive.google.com/file/d/1aMYxG5YJUgMHpbWmHpzL4hRfmP26UQqh) | [step 5](https://drive.google.com/file/d/17kQt2r78ReWle4PhEIOXG7w7BFdezGM1) | [summary](https://drive.google.com/file/d/1639IFUG7poaDE2uB6aISUOi65ooBwCIg) |
| ------------------ | --------------------------------------------------------------------------- | --------------------------------------------------------------------------- | --------------------------------------------------------------------------- | --------------------------------------------------------------------------- | --------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |

### Halle University

> Edit 2024-05-31: These links are now dead. The data is still available on Materials Cloud and Google Drive.

The [paper itself][wbm paper] links to a discontinued site `https://tddft.org/bmg/data.php` which listed download URLs for CIF files and the `ComputedStructureEntries` (CSEs) of steps 1-3:

| Halle University links | [step 1 CSEs](https://tddft.org/bmg/files/data/substitutions_000.json.bz2) | [step 2 CSEs](https://tddft.org/bmg/files/data/substitutions_001.json.bz2) | [step 3 CSEs](https://tddft.org/bmg/files/data/substitutions_002.json.bz2) | [CIF files](https://tddft.org/bmg/files/data/similarity-cifs.tar.gz) |
| ---------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------- |

### Materials Cloud

materialscloud:2021.68 includes a readme describing the dataset and the summary CSV columns, plus a Python script for loading the data.

| [Materials Cloud](https://archive.materialscloud.org/record/2021.68) | [step 1](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/step_1.json.bz2/content) | [step 2](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/step_2.json.bz2/content) | [step 3](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/step_3.json.bz2/content) | [step 4](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/step_4.json.bz2/content) | [step 5](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/step_5.json.bz2/content) | [summary](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/summary.txt.bz2/content) | [readme](https://archive.materialscloud.org/api/records/qp8qv-gry23/files/README.txt/content) |
| -------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |

[wbm paper]: https://nature.com/articles/s41524-020-00481-6

## 🎯 &thinsp; Target Distribution

Below is a histogram of hull distances of all WBM test set structures with respect to the MP convex hull, the target distribution the models are asked to predict. The dummy MAE of always predicting the test set mean is **0.12 eV/atom**.

The number of stable materials (according to the MP convex hull spanned by the training data the models have access to) is **42.7k** out of **257k**, a dummy stability hit rate of **16.6%**.

> Note: [According to the WBM authors](https://nature.com/articles/s41524-020-00481-6#Sec2), the stability rate w.r.t. the more complete hull constructed from the combined train and test set (MP + WBM) for the first 3 rounds of elemental substitution is 18,479 out of 189,981 crystals ($\approx$ 9.7%), so many WBM structures are stable w.r.t. MP but not w.r.t. each other. One could build an active learning benchmark from this dataset where a well-calibrated model should initially predict a ~16.6% stability rate but converge towards 9.7% as it is retrained on more and more chunks of the WBM test set.

{@render hist_wbm_hull_dist?.()}

## 🧪 &thinsp; Chemical Diversity

The WBM test set and even more so the MP training set are heavily oxide dominated. The WBM test set is about 66% larger than the MP training set and more chemically diverse, with a higher fraction of transition metals, post-transition metals and metalloids. We picked such a large, diverse test set for future-proofing: ideally, it will remain a challenging materials discovery test bed even for large foundational ML models.

### WBM

The WBM test set consists of 256,963 WBM `ComputedStructureEntries`

{@render wbm_elements_heatmap?.()}

### MP

The MP training set consists of 154,718 `ComputedStructureEntries`

{@render mp_elements_heatmap?.()}

### MPtrj

{@render mp_trj_elements_heatmap?.()}

## 📊 &thinsp; Symmetry Statistics

These sunburst diagrams show the spacegroup distribution of MP (left) and WBM (right). Both cover all 7 crystal systems well, except triclinic crystals which are just 1% of WBM but well represented in MP (15%). The 3 largest systems are monoclinic, orthorhombic and triclinic in MP vs orthorhombic, tetragonal and cubic in WBM. WBM structures thus have overall higher symmetry, which can benefit some models more than others. Wrenformer in particular uses symmetries as a coarse-grained description of the underlying structure, so its representation essentially degrades to composition only on symmetry-less P1 structures. Given this spacegroup distribution, it should fare well on the WBM test set. That it is still outperformed by all interatomic potentials and some single-shot GNNs indicates the underlying methodology is unable to compete.

{@render spacegroup_sunbursts?.()}

## Downloading WBM data

The `DataFiles` class in `matbench_discovery.data` automatically downloads and locally caches data files needed for model training and testing. For example, to download and access the WBM summary data:

```py
from matbench_discovery.data import DataFiles, ase_atoms_from_zip
import pandas as pd

df_wbm = pd.read_csv(DataFiles.wbm_summary.path)
# confirm test set size
assert df_wbm.shape == (256_963, 18)

# WBM initial structures in pymatgen JSON format
df_init_structs = pd.read_json(DataFiles.wbm_initial_structures.path, lines=True)
# Columns: material_id, formula_from_cse, initial_structure
# WBM initial structures as ASE Atoms
wbm_init_atoms = ase_atoms_from_zip(DataFiles.wbm_initial_atoms.path)
assert len(wbm_init_atoms) == 256_963
```

`"wbm-summary"` columns:

1. **`material_id`**: Structure identifier encoding the substitution step and index, as described above.
1. **`formula`**: A compound's unreduced alphabetical formula
1. **`n_sites`**: Number of sites in the structure's unit cell
1. **`volume`**: Relaxed structure volume in cubic Angstrom
1. **`uncorrected_energy`**: Raw VASP-computed energy
1. **`e_form_per_atom_wbm`**: Original formation energy per atom from [WBM paper]
1. **`e_above_hull_wbm`**: Original energy above the convex hull (eV/atom) from [WBM paper]
1. **`protostructure_spglib`**: Aflow label strings built from spacegroup and Wyckoff positions of the DFT-relaxed structure as computed by [spglib](https://spglib.readthedocs.io/en/stable/api/autodoc/spglib.spg.html#spglib.spg.get_symmetry_dataset).
1. **`protostructure_spglib_initial_structure`**: Same as `protostructure_spglib` but computed from the initial structure.
1. **`bandgap_pbe`**: PBE-level DFT band gap from [WBM paper]
1. **`uncorrected_energy_from_cse`**: Uncorrected DFT energy stored in `ComputedStructureEntries`. Should equal `uncorrected_energy`. In 2 cases (`wbm-2-3218`, `wbm-1-56320`) the absolute difference between the summary file and the computed structure entries exceeds 0.1 eV (0.14 eV each), too large for the summary's 3-decimal rounding. Their cause is unknown, perhaps data corruption unnoticed by the WBM authors.
1. **`e_form_per_atom_uncorrected`**: Uncorrected DFT formation energy per atom in eV/atom.
1. **`e_form_per_atom_mp2020_corrected`**: Matbench Discovery's ground truth for the formation energy: the [MP2020 energy corrections][`MaterialsProject2020Compatibility`] (current `pymatgen` defaults, see the correction note above) applied to `e_form_per_atom_uncorrected`.
1. **`e_correction_per_atom_mp2020`**: [`MaterialsProject2020Compatibility`] energy corrections in eV/atom.
1. **`e_correction_per_atom_mp_legacy`**: Legacy [`MaterialsProjectCompatibility`] energy corrections in eV/atom. Having both old and new corrections allows updating predictions from older models like MEGNet trained on MP formation energies with the old correction scheme. Carried over unchanged from the 2023-12-13 summary.
1. **`e_above_hull_mp2020_corrected_ppd_mp`**: Energy above hull in eV/atom after applying the MP2020 correction scheme, w.r.t. the convex hull spanned by all 154,718 Materials Project `ComputedStructureEntries`. Matbench Discovery takes these as ground truth for material stability. Any value above 0 is assumed to be an unstable/metastable material.
1. **`site_stats_fingerprint_init_final_norm_diff`**: Norm of the difference between the initial and final site fingerprints, a volume-independent measure of how much the structure changed during DFT relaxation. Uses the `matminer` [`SiteStatsFingerprint`](https://github.com/hackingmaterials/matminer/blob/33bf1120/matminer/featurizers/structure/sites.py#L21-L33) (v0.8.0).
1. **`unique_prototype`**: Flag selecting the unique-prototype subset described above.

[`MaterialsProject2020Compatibility`]: https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/compatibility.py#L823
[`MaterialsProjectCompatibility`]: https://github.com/materialsproject/pymatgen/blob/02a4ca8aa/pymatgen/entries/compatibility.py#L766
