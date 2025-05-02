<script lang="ts">
  import StructRmsdCdfModels from '$figs/struct-rmsd-cdf-models.svelte'
  import SymOpsDiffBar from '$figs/sym-ops-diff-bar.svelte'
  import { GeoOptMetricsTable, MODELS } from '$lib'
  import { pretty_num } from 'elementari'
  import type { SvelteComponent } from 'svelte'
  import GeoOptReadme from './geo-opt-readme.md'

  const plots = import.meta.glob(`$figs/spg-sankey-*symprec=1e-05.svelte`, {
    eager: true,
    import: `default`,
  }) as Record<string, typeof SvelteComponent>

  const n_min_relaxed_structures: number = MODELS.reduce((acc, model) => {
    const geo_opt = model.metrics?.geo_opt
    if (typeof geo_opt === `string` || !geo_opt) return acc
    return Math.min(acc, geo_opt?.[`symprec=1e-2`]?.n_structures ?? Infinity)
  }, Infinity)
</script>

<GeoOptReadme>
  {#snippet geo_opt_metrics_table()}
    <GeoOptMetricsTable show_non_compliant />
  {/snippet}
  {#snippet min_relaxed_structures()}
    <span>{pretty_num(n_min_relaxed_structures)}</span>
  {/snippet}
  {#snippet struct_rmsd_cdf_models()}
    <StructRmsdCdfModels />
  {/snippet}
  {#snippet sym_ops_diff_bar()}
    <SymOpsDiffBar />
  {/snippet}
</GeoOptReadme>

<ul>
  {#each Object.entries(plots) as [name, Plot], idx (name + idx)}
    <Plot {name} style="width: 100%; place-self: center;" />
  {/each}
</ul>

<style>
  ul {
    padding: 0;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
    gap: 3em 2em;
  }
</style>
