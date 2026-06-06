<script lang="ts">
  import ref_energies from '$figs/mp-elemental-ref-energies.json.gz'
  import { ScatterPlot } from 'matterviz/plot'
</script>

## MP Elemental Reference Energies

<ScatterPlot
series={[{
...ref_energies,
markers: `line+points`,
line_style: { stroke: `#636efa` },
point_style: { fill: `#636efa` },
}]}
x_axis={{ label: `Atomic number`, range: [0, null] }}
y_axis={{ label: `Energy (eV/atom)` }}
legend={null}
style="height: 420px"
/>

<!-- > @label:fig:mp-elemental-reference-energies -->

The above MP elemental reference energies are the lowest energies found for unary structures for that element in MP.
These values are used to calculate the formation energy target values of WBM materials by subtracting from each WBM material's DFT energy the sum of the elemental reference energies multiplied by their stoichiometric coefficients.<br>
These database entries were [queried on 2023-02-07](https://github.com/janosh/matbench-discovery/blob/-/data/mp/2023-02-07-mp-elemental-reference-entries.json.gz).
