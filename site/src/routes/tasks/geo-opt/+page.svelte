<script lang="ts">
  import { browser } from '$app/environment'
  import StructRmsdCdfModels from '$figs/struct-rmsd-cdf-models.svelte'
  import SymOpsDiffBar from '$figs/sym-ops-diff-bar.svelte'
  import { GeoOptMetricsTable } from '$lib'
  import type { SvelteComponent } from 'svelte'
  import GeoOptReadme from './geo-opt-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>
</script>

<GeoOptReadme>
  <GeoOptMetricsTable show_non_compliant slot="geo-opt-metrics-table" />
  <StructRmsdCdfModels slot="struct-rmsd-cdf-models" />
  <SymOpsDiffBar slot="sym-ops-diff-bar" />
</GeoOptReadme>

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
