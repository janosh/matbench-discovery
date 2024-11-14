<script lang="ts">
  import { browser } from '$app/environment'
  import StructRmsdCdfModels from '$figs/struct-rmsd-cdf-models.svelte'
  import SymOpsDiffBar from '$figs/sym-ops-diff-bar.svelte'
  import { GeoOptMetricsTable } from '$lib'
  import type { SvelteComponent } from 'svelte'
  import SankeyReadme from './sankey-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>
</script>

<h1>MLFF Geometry Optimization Analysis</h1>

All plots/metrics below evaluate the quality of MLFF relaxations for the 257k crystal
structures in the WBM test set. Not all models were able to relax all structures
(user/cluster error may explain some failures) but every model was evaluated on at least
240k relaxations.

<figure>
  <GeoOptMetricsTable show_non_compliant />
  <figcaption>
    σ<sub>match</sub> / σ<sub>dec</sub> / σ<sub>inc</sub> denote the fraction of
    structures that retain, increase, or decrease the symmetry of the DFT-relaxed
    structure during MLFF relaxation. The match criterion is for the ML ground state to
    have identical spacegroup as DFT (according to spglib with default precision
    settings). For σ<sub>dec</sub> / σ<sub>inc</sub>, ML relaxation increased / decreased
    the set of symmetry operations on a structure.
  </figcaption>
</figure>

<hr />

{#if browser}
  <figure>
    <StructRmsdCdfModels />
    <figcaption>
      Cumulative distribution of RMSD between ML and DFT-relaxed structures.
    </figcaption>
  </figure>
{/if}

<hr />

{#if browser}
  <figure>
    <SymOpsDiffBar />
    <figcaption>
      Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models
      are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub
        >sym ops,ML</sub
      >
      - N<sub>sym ops,DFT</sub>
      distribution.
    </figcaption>
  </figure>
{/if}

<hr />

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
