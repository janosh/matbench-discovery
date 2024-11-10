<script lang="ts">
  import { browser } from '$app/environment'
  import { GeoOptMetricsTable } from '$lib'
  import type { SvelteComponent } from 'svelte'
  import SankeyReadme from './sankey-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>
</script>

<GeoOptMetricsTable show_non_compliant />

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
