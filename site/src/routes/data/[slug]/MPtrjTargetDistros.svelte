<script lang="ts">
  import type { Component } from 'svelte'
  import MPtrjElemCountsPtable from './MPtrjElemCountsPtable.svelte'
  import MPtrjTargetCounts from './mptrj-target-counts.md'

  const plots = import.meta.glob<Component>(
    `$figs/mp-trj-*.svelte`,
    { eager: true, import: 'default' },
  )

  const title_map: Record<string, string> = {
    'e-form': `Formation Energy`,
    forces: `Forces`,
    stresses: `Stresses`,
    magmoms: `Magnetic Moments`,
  }
</script>

<MPtrjTargetCounts />

<ul>
  {#each Object.entries(plots) as [name, Plot] (name)}
    {@const title = name.split(`mp-trj-`)[1].split(`-hist.svelte`)[0]}
    <div>
      <h3>{title_map[title] ?? title}</h3>
      <Plot {title} style="width: 100%; max-width: 700px; height: 300px" />
    </div>
  {/each}
</ul>

<h2>Elemental Prevalence</h2>

<MPtrjElemCountsPtable />

<style>
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    gap: 1em;
  }
  ul h3 {
    text-transform: capitalize;
    text-align: center;
    margin: 1em auto 0;
  }
</style>
