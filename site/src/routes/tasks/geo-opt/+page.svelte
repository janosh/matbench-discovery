<script lang="ts">
  import { browser } from '$app/environment'
  import StructRmsdCdfModels from '$figs/struct-rmsd-cdf-models.svelte'
  import { GeoOptMetricsTable } from '$lib'
  import type { SvelteComponent } from 'svelte'
  import SankeyReadme from './sankey-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>
</script>

<figure>
  <GeoOptMetricsTable show_non_compliant />
  <figcaption>
    Match / Increase / Decrease count structures that retain, increase, or decrease the
    symmetry of the DFT-relaxed structure. The match criterion is for the ML ground state
    to have identical spacegroup as DFT (according to spglib with default precision
    settings). Increase / decrease mean the set of symmetry operations on the structure
    grew / shrank.
  </figcaption>
</figure>

<figure>
  <StructRmsdCdfModels />
  <figcaption>
    Cumulative distribution of RMSD comparing ML vs DFT-relaxed structures for each model.
  </figcaption>
</figure>

<SankeyReadme />

{#if browser}
  <ul>
    {#each Object.entries(plots) as [name, Plot]}
      <Plot {name} style="width: 100%; place-self: center;" />
    {/each}
  </ul>
{/if}

<style>
  ul {
    padding: 0;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
    gap: 3em 2em;
  }
</style>
