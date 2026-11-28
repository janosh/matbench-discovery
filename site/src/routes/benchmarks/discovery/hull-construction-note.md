### Convex Hull Construction in Matbench Discovery

In Matbench Discovery, **the convex hull is always constructed from DFT reference energies**, not from the ML model's predicted energies. This methodological choice differs from some other benchmarking approaches and is key to correctly interpreting the energy metrics.

#### What This Means

- **DFT-based hull:** To calculate a material's distance to the convex hull (E<sub>hull dist</sub>), we compare the model's predicted formation energy against the DFT convex hull built from Materials Project reference structures.
- **Fixed reference:** The hull does not depend on the model's predictions. All models are evaluated against the same DFT reference hull.
- **Discovery criterion:** A material counts as a "discovery" if the model correctly predicts it to be on or below the convex hull spanned by all known DFT-computed competing phases in its chemical system in Materials Project. The MP reference entries were queried on 2023-02-07, database release [v2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28).

#### Why This Matters

1. **Formation energy MAE = Hull distance MAE:** The model's prediction and the DFT reference share the same energy scale (relative to the same elemental references), and the fixed hull shifts predicted and true hull distances by the same per-material offset. The hull distance error therefore equals the formation energy error for every material.

2. **Systematic errors are not canceled:** Systematic model errors (e.g. consistently over- or underpredicting certain elements) appear in both the formation energy and hull distance metrics. The model cannot "correct" for them by having them affect the test structures and the reference hull equally.
   - **Advantage:** Tests absolute accuracy of model predictions against ground truth.
   - **Use case:** Pre-screening candidates for DFT calculations and evaluating a model's ability to identify materials below the DFT reference hull.

3. **Different from some literature:** Papers like [Bartel et al., npj Computational Materials 6, 97 (2020)](https://www.nature.com/articles/s41524-020-00362-y) construct hulls from model predictions, allowing systematic model errors to partially cancel. There, the hull distance MAE can differ from the formation energy MAE.
   - **Advantage:** Systematic model errors can partially cancel.
   - **Use case:** Using the model as a complete replacement for DFT.

This subtle distinction matters for interpreting model performance and comparing benchmarks fairly.
