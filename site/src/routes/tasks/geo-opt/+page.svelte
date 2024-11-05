<script lang="ts">
  import { browser } from '$app/environment'
  import type { SvelteComponent } from 'svelte'
  import SankeyReadme from './sankey-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>
</script>

<SankeyReadme />

{#if browser}
  <ul>
    {#each Object.entries(plots) as [name, Plot]}
      <Plot {name} style="width: 100%; max-width: 700px; max-height: 400px;" />
    {/each}
  </ul>
{/if}

<style>
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    gap: 7em;
  }
</style>
