# WBM Dataset

Source: [Predicting stable crystalline compounds using chemical similarity](https://nature.com/articles/s41524-020-00481-6) (2021)


## `wbm-summary.csv`

Load with

```py
df_wbm_summary = pd.read_csv(  # download wbm-steps-summary.csv (23.31 MB)
    "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")
```

## Comprehensive Link Collection for WBM dataset

Links to WBM data files have proliferated. This is an attempt to keep track of all of them.

Initial structures were sent as Google Drive links via email by Hai-Chen Wang on 2021-09-01.

step 1: https://drive.google.com/file/d/1ZUgtYwrfZn_P8bULWRtTXepyAxHVxS5C
step 2: https://drive.google.com/file/d/1-3uu2AcARJxH7GReteGVASZTuttFGiW_
step 3: https://drive.google.com/file/d/1hc5BvDiFfTu_tc5F8m7ONSw2OgL9vN6o
step 4: https://drive.google.com/file/d/1aMYxG5YJUgMHpbWmHpzL4hRfmP26UQqh
step 5: https://drive.google.com/file/d/17kQt2r78ReWle4PhEIOXG7w7BFdezGM1
summary: https://drive.google.com/file/d/1639IFUG7poaDE2uB6aISUOi65ooBwCIg

The `ComputedStructureEntries` for steps 1-3 were also linked from the Nature paper:

Index page: https://tddft.org/bmg/data.php
step 1 CSEs: https://tddft.org/bmg/files/data/substitutions_000.json.bz2
step 2 CSEs: https://tddft.org/bmg/files/data/substitutions_001.json.bz2
step 3 CSEs: https://tddft.org/bmg/files/data/substitutions_002.json.bz2
CIF files: https://tddft.org/bmg/files/data/similarity-cifs.tar.gz

Materials Cloud archive: https://archive.materialscloud.org/record/2021.68
File URLs:
readme: https://archive.materialscloud.org/record/file?record_id=840&filename=README.txt
summary: https://archive.materialscloud.org/record/file?record_id=840&filename=summary.txt.bz2
step 1: https://archive.materialscloud.org/record/file?record_id=840&filename=step_1.json.bz2 etc.
