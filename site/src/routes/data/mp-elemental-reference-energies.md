<script>
  import MPRefEnergies from '$figs/mp-elemental-ref-energies.svelte'
</script>

## MP Elemental Reference Energies

<MPRefEnergies />

<!-- > @label:fig:mp-elemental-reference-energies -->

The above MP elemental reference energies are the lowest energies found for unary structures for that element in MP.
These values are used to calculate the formation energy target values of WBM materials by subtracting from each WBM material's DFT energy the sum of the elemental reference energies multiplied by their stoichiometric coefficients.<br>
These database entries were [queried on 2023-02-07](https://github.com/janosh/matbench-discovery/blob/-/data/mp/2023-02-07-mp-elemental-reference-entries.json.gz). Marker size indicates the number of atoms in the reference structure. Hover for details like material ID and energy per atom.
