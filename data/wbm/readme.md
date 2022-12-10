# WBM Dataset

The **WBM dataset** was published in [Predicting stable crystalline compounds using chemical similarity][wbm paper] (Nature Computational Materials, Jan 2021, [doi:10.1038/s41524-020-00481-6](http://doi.org/10.1038/s41524-020-00481-6)). The authors generated 257,487 structures through single-element substitutions on Materials Project (MP) source structures. The replacement element was chosen based on chemical similarity determined by a matrix data-mined from the [Inorganic Crystal Structure Database (ICSD)](https://icsd.products.fiz-karlsruhe.de).

The resulting novel structures were relaxed using MP-compatible VASP inputs (i.e. using `pymatgen`'s `MPRelaxSet`) and identical POTCARs in an attempt to create a database of Materials Project compatible novel crystals. Any degrade in model performance from training to test set should therefore largely be a result of extrapolation error rather than covariate shift in the underlying data.

The authors performed 5 rounds of elemental substitution in total, each time relaxing all generated structures and adding those found to lie on the convex hull back to the source pool. In total, ~20k or close to 10% were found to lie on the Materials Project convex hull.

Since repeated substitutions should - on average - increase chemical dissimilarity, the 5 iterations of this data-generation process are a unique and compelling feature as it allows out-of distribution testing. We can check how model performance degrades when asked to predict on structures increasingly more dissimilar from the training set (which is restricted to the MP 2022 database release (or earlier) for all models in this benchmark).

## About the IDs

As you may have guessed, the first integer in each material ID following the prefix `wbm-` ranges from 1 to 5 and indicates the substitution iteration count. Each iteration has varying numbers of materials counted by the 2nd integer. Note the 2nd integer is not strictly consecutive. A small number of materials (~0.2%) were removed by the data processing steps detailed below. Don't be surprised to find an ID like `wbm-3-70804` followed by

## Data processing steps

The full set of processing steps used to curate the WBM test set from the raw data files (downloaded from the URLs listed below) can be found in [`data/wbm/fetch_process_wbm_dataset.py`](https://github.com/janosh/matbench-discovery/blob/site/data/wbm/fetch_process_wbm_dataset.py). Processing involved

- re-formatting material IDs
- correctly aligning initial structures to DFT-relaxed `ComputedStructureEntries`
- remove 6 pathological structures (with 0 volume)
- remove formation energy outliers below -5 and above 5 eV/atom (removed 502 and 22 crystals respectively out of 257,487 total, including an anomaly of 500 structures at exactly -10 eV/atom)
  <!-- ![WBM formation energy histogram indicating outlier cutoffs](2022-12-07-hist-e-form-per-atom.png) -->
- apply the [`MaterialsProject2020Compatibility`](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility) energy correction scheme to the formation energies
- compute energy to the convex hull constructed from all MP `ComputedStructureEntries` queried on 2022-09-16 ([database release 2021.05.13](https://docs.materialsproject.org/changes/database-versions#v2021.05.13))

The number of materials in each step before and after processing are:

| step | 1      | 2      | 3      | 4      | 5      | total   |
| ---- | ------ | ------ | ------ | ------ | ------ | ------- |
| pre  | 61,848 | 52,800 | 79,205 | 40,328 | 23,308 | 257,487 |
| post | 61,466 | 52,755 | 79,160 | 40,314 | 23,268 | 256,963 |

Invoking that script with `python fetch_process_wbm_dataset.py` will auto-download and regenerate the WBM test set files from scratch. If you find any questionable in the released test set or inconsistencies between the files on GitHub vs the output of that script, please [raise an issue](https://github.com/janosh/matbench-discovery/issues).

## Links to WBM data files

Links to WBM data files have proliferated. This is an attempt to keep track of all of them.

Initial structures were sent as Google Drive links via email by Hai-Chen Wang on 2021-09-01.

step 1: <https://drive.google.com/file/d/1ZUgtYwrfZn_P8bULWRtTXepyAxHVxS5C>
step 2: <https://drive.google.com/file/d/1-3uu2AcARJxH7GReteGVASZTuttFGiW_>
step 3: <https://drive.google.com/file/d/1hc5BvDiFfTu_tc5F8m7ONSw2OgL9vN6o>
step 4: <https://drive.google.com/file/d/1aMYxG5YJUgMHpbWmHpzL4hRfmP26UQqh>
step 5: <https://drive.google.com/file/d/17kQt2r78ReWle4PhEIOXG7w7BFdezGM1>
summary: <https://drive.google.com/file/d/1639IFUG7poaDE2uB6aISUOi65ooBwCIg>

The `ComputedStructureEntries` for steps 1-3 were also linked from the [WBM Nature paper][wbm paper]:

Index page: <https://tddft.org/bmg/data.php>
step 1 CSEs: <https://tddft.org/bmg/files/data/substitutions_000.json.bz2>
step 2 CSEs: <https://tddft.org/bmg/files/data/substitutions_001.json.bz2>
step 3 CSEs: <https://tddft.org/bmg/files/data/substitutions_002.json.bz2>
CIF files: <https://tddft.org/bmg/files/data/similarity-cifs.tar.gz>

Materials Cloud archive: <https://archive.materialscloud.org/record/2021.68>
File URLs:

- readme: <https://archive.materialscloud.org/record/file?record_id=840&filename=README.txt>
- summary: <https://archive.materialscloud.org/record/file?record_id=840&filename=summary.txt.bz2>
- step 1: <https://archive.materialscloud.org/record/file?record_id=840&filename=step_1.json.bz2>
- step 2: <https://archive.materialscloud.org/record/file?record_id=840&filename=step_2.json.bz2>
- step 3: <https://archive.materialscloud.org/record/file?record_id=840&filename=step_3.json.bz2>
- step 4: <https://archive.materialscloud.org/record/file?record_id=840&filename=step_4.json.bz2>
- step 5: <https://archive.materialscloud.org/record/file?record_id=840&filename=step_5.json.bz2>

[wbm paper]: https://nature.com/articles/s41524-020-00481-6
