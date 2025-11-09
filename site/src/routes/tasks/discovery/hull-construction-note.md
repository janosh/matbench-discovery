## Convex Hull Construction in Matbench Discovery

In Matbench Discovery, **the convex hull is always constructed from DFT reference energies**, not from the ML model's predicted energies. This is an important methodological choice that differs from some other benchmarking approaches and has several implications. Understanding how the convex hull is constructed is important for correctly interpreting the energy metrics in Matbench Discovery.

### What This Means

- **DFT-based hull:** When we calculate the distance to the convex hull (E<sub>hull dist</sub>) for a material, we compare the model's predicted formation energy against the DFT-computed convex hull built from Materials Project reference structures.
- **Fixed reference:** The hull does not change based on the model's predictions. All models are evaluated against the same DFT reference hull.
- **Discovery criterion:** A material is counted as a "discovery" if the model correctly predicts it to be lower in energy than all known DFT-computed competing phases with the same (reduced) composition in Materials Project. The reference data was pulled on 2023-03-16 (14 GB), database release [v2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28).

### Why This Matters

This approach means that:

1. **Formation energy MAE = Hull distance MAE:** Because both the model's prediction and the DFT reference are measured on the same energy scale (relative to the same elemental references), the error in formation energy prediction directly equals the error in hull distance prediction. This is a consequence of linear transformations leaving the MAE metric invariant.

2. **Systematic errors are not canceled:** If a model has systematic errors (e.g., consistently over- or underpredicting certain elements), these errors will appear in both the formation energy and hull distance metrics. The model cannot "correct" for its own systematic errors by having them affect both the test structures and the reference hull equally.

   - **Advantage:** Tests absolute accuracy of model predictions against ground truth.
   - **Use case:** Pre-screening candidates for DFT calculations and evaluating a model's ability to identify materials below the DFT reference hull.

3. **Different from some literature:** Papers like [Nature Communications 11:3793 (2020)](https://www.nature.com/articles/s41524-020-00362-y) construct hulls from model predictions, allowing systematic model errors to partially cancel. In that approach, the hull distance MAE can differ from the formation energy MAE.

   - **Advantage:** Systematic model errors can partially cancel.
   - **Use case:** Using the model as a complete replacement for DFT.

This distinction is subtle but important for correctly interpreting model performance and making fair comparisons between different benchmarks.
