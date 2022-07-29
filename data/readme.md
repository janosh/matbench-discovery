# Data file descriptions

## `wbm-cses-and-initial-structures.json.gz`

257486 rows × 3 columns:

- cse: pymatgen `ComputedStructureEntries`
- initial_structure: pymatgen `Structure` objects
- bandgap: `float`

## `wbm-steps-summary.csv`

Generated in `dielectric_frontier` project.

Columns: `material_id, formula, n_sites, volume, energy, e_form, e_hull, e_above_hull, bandgap_pbe`

Load with

```py
df_summary = pd.read_csv(f"{ROOT}/data_others/wbm/data/steps-summary.csv", comment="#").set_index("material_id")
```
