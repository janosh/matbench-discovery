# Data file descriptions

## `wbm-cses-and-initial-structures.json.gz`

Source: [Predicting stable crystalline compounds using chemical similarity](https://nature.com/articles/s41524-020-00481-6) (2021)

257486 rows × 3 columns:

- cse: pymatgen `ComputedStructureEntries`
- initial_structure: pymatgen `Structure` objects
- bandgap: `float`

## `wbm-steps-summary.csv`

Generated in `dielectric_frontier` project.

Columns: `material_id, formula, n_sites, volume, energy, e_form, e_hull, e_above_hull, bandgap_pbe`

Load with

```py
df_wbm_summary = pd.read_csv(  # download wbm-steps-summary.csv (23.31 MB)
    "https://figshare.com/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")
```
