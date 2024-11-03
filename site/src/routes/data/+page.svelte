<script lang="ts">
  import { browser } from '$app/environment'
  import FormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import HistWbmHullDist from '$figs/hist-wbm-hull-dist.svelte'
  import MPtrjNSitesHist from '$figs/mp-trj-n-sites-hist.svelte'
  import MPvsMPtrjVsWBMArityHist from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.svelte'
  import SpacegroupSunburstMp from '$figs/spacegroup-sunburst-mp.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import { PtableHeatmap } from '$lib'
  import DataReadme from '$root/data/wbm/readme.md'
  import { ColorScaleSelect } from 'elementari'
  import Select from 'svelte-multiselect'
  import type { Snapshot } from './$types'
  import MpElementalReferenceEnergies from './mp-elemental-reference-energies.md'
  import MPtrjElemCountsPtable from './mptrj/MPtrjElemCountsPtable.svelte'

  const elem_counts = import.meta.glob(
    `./*-element-counts-by-{occurrence,composition}*.json`,
    {
      eager: true,
      import: `default`,
    },
  )

  let log = false // log color scale
  let color_scale = [`Viridis`]
  const count_modes = [`occurrence`, `composition`]
  let count_mode = count_modes[0]

  $: mp_elem_counts = elem_counts[`./mp-element-counts-by-${count_mode}.json`]
  $: if (!mp_elem_counts) throw `No MP data for count mode ${count_mode}!`
  $: mp_trj_elem_counts = elem_counts[`./mp-trj-element-counts-by-${count_mode}.json`]
  $: if (!mp_trj_elem_counts) throw `No MPtrj data for count mode ${count_mode}!`
  $: wbm_elem_counts = elem_counts[`./wbm-element-counts-by-${count_mode}.json`]
  $: if (!wbm_elem_counts) throw `No WBM data for count mode ${count_mode}!`

  export const snapshot: Snapshot = {
    capture: () => ({ color_scale, log, count_mode }),
    restore: (values) => ({ color_scale, log, count_mode } = values),
  }
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if browser}
      <FormEnergyHist />
    {/if}
  </svelte:fragment>
  <svelte:fragment slot="wbm-elements-heatmap">
    <label
      for="count-mode"
      style="display: inline-block; transform: translate(10cqw, 5ex);">Count Mode</label
    >
    <Select
      id="count-mode"
      selected={[count_mode]}
      bind:value={count_mode}
      options={count_modes}
      minSelect={1}
      maxSelect={1}
    />
    <ColorScaleSelect bind:selected={color_scale} />
    <p style="text-align: center; margin: 2em; font-size: smaller;">
      The difference between count modes is best explained by example.
      <code>occurrence</code> mode maps Fe<sub>2</sub>O<sub>3</sub> to {`{Fe: 1, O: 1}`},
      <code>composition</code>
      mode maps it to {`{Fe: 2, O: 3}`}.
    </p>
    <PtableHeatmap
      heatmap_values={wbm_elem_counts}
      color_scale={color_scale[0]}
      color_bar_props={{ label: `WBM element counts by ${count_mode}` }}
      {log}
    />
  </svelte:fragment>

  <svelte:fragment slot="mp-elements-heatmap">
    <PtableHeatmap
      heatmap_values={mp_elem_counts}
      color_scale={color_scale[0]}
      color_bar_props={{ label: `MP element counts by ${count_mode}` }}
      {log}
    />
  </svelte:fragment>

  <svelte:fragment slot="mp-trj-elements-heatmap">
    <MPtrjElemCountsPtable {count_mode} {log} color_scale={color_scale[0]} />
  </svelte:fragment>

  <svelte:fragment slot="hist-wbm-hull-dist">
    {#if browser}
      <HistWbmHullDist />
    {/if}
  </svelte:fragment>

  <div
    style="display: flex; gap: 1em; justify-content: space-around; flex-wrap: wrap;"
    slot="spacegroup-sunbursts"
  >
    {#if browser}
      <SpacegroupSunburstMp />
      <SpacegroupSunburstWbm />
    {/if}
  </div>
</DataReadme>

<MpElementalReferenceEnergies />

{#if browser}
  <MPvsMPtrjVsWBMArityHist style="margin: auto; max-width: 60cqw; padding-right: 2em;" />
{/if}
<p>
  Distribution of unique elements per structure in MP, MPtrj and WBM. The bar heights are
  normalized by the total number of structures in each data set. WBM is dominated by
  ternary phases making up 74% of the data set followed by about 13% each of binaries and
  quaternaries. MP has a more even distribution, in particular with more than double the
  relative share of quaternary phases and a significant number of quinternaries which are
  almost absent from WBM. Not shown in this plot for visual clarity are 3% of MP
  structures containing more than 5 elements (up to 9). We also include MPtrj in this plot
  to show a slight drop in relative abundance of quinternary and higher phases vs MP
  ground states.
</p>

<MPtrjNSitesHist style="margin: auto; max-width: 80cqw; padding-right: 2em;" />
<p>
  Histogram of number of atoms per structure. The inset shows the same distribution
  log-scaled to visualize the tail of large structures. The green cumulative line in the
  inset shows that 82% have less than 50 sites and 97% of structures in MPtrj have less
  than 100 atoms.
</p>

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
    justify-items: center;
  }
</style>
