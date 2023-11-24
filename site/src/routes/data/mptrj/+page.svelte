<script lang="ts">
  import { browser } from '$app/environment'
  import type { SvelteComponent } from 'svelte'
  import MPtrjElemCountsPtable from './MPtrjElemCountsPtable.svelte'
  import Readme from './readme.md'

  const plots = import.meta.glob(`$figs/mp-trj-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>

  const title_map: Record<string, string> = {
    'e-form': `Formation Energy`,
    forces: `Forces`,
    stresses: `Stresses`,
    magmoms: `Magnetic Moments`,
  }
</script>

<Readme />

{#if browser}
  <ul>
    {#each Object.entries(plots) as [name, Plot]}
      {@const title = name.split(`mp-trj-`)[1].split(`-hist.svelte`)[0]}
      <h2>{title_map[title] ?? title}</h2>
      <Plot {title} style="width: 100%; max-width: 700px; max-height: 400px;" />
    {/each}
  </ul>
{/if}

<h2>Elemental Prevalence</h2>

<MPtrjElemCountsPtable />

<style>
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    gap: 1em;
  }
  ul h2 {
    text-transform: capitalize;
    margin: 1em auto 0;
  }
</style>
