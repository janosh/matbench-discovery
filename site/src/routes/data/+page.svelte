<script lang="ts">
  import { browser } from '$app/environment'
  import FormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import HistWbmHullDist from '$figs/hist-wbm-hull-dist.svelte'
  import SpacegroupSunburstMp from '$figs/spacegroup-sunburst-mp.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import { PtableInset } from '$lib'
  import DataReadme from '$root/data/wbm/readme.md'
  import type { ChemicalElement } from 'elementari'
  import { ColorBar, ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import Select from 'svelte-multiselect'
  import { Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'
  import MpElementalReferenceEnergies from './mp-elemental-reference-energies.md'

  const elem_counts = import.meta.glob(`./*-element-counts-{occu,comp}*.json`, {
    eager: true,
    import: `default`,
  })

  let log = false // log color scale
  let color_scale = [`Inferno`]
  let active_mp_elem: ChemicalElement
  let active_wbm_elem: ChemicalElement
  const count_mode_ops = [`occurrence`, `composition`]
  let count_mode = count_mode_ops[0]

  $: mp_elem_counts = elem_counts[`./mp-element-counts-${count_mode}.json`]
  $: wbm_elem_counts = elem_counts[`./wbm-element-counts-${count_mode}.json`]

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
    <PeriodicTable
      heatmap_values={wbm_elem_counts}
      color_scale={color_scale[0]}
      {log}
      bind:active_element={active_wbm_elem}
    >
      <TableInset slot="inset">
        <label for="log">Log color scale<Toggle id="log" bind:checked={log} /></label>
        <PtableInset element={active_wbm_elem} elem_counts={wbm_elem_counts} />
        <ColorBar
          text_side="top"
          color_scale={color_scale[0]}
          tick_labels={5}
          precision={3}
          range={[0, Math.max(...Object.values(wbm_elem_counts))]}
          style="width: 85%; margin: 0 2em 2em;"
        />
      </TableInset>
    </PeriodicTable>
    <label
      for="count-mode"
      style="display: inline-block; transform: translate(10cqw, 5ex);">Count Mode</label
    >
    <Select
      id="count-mode"
      selected={[count_mode]}
      bind:value={count_mode}
      options={count_mode_ops}
      minSelect={1}
      maxSelect={1}
    />
    <ColorScaleSelect bind:selected={color_scale} />
  </svelte:fragment>

  <svelte:fragment slot="mp-elements-heatmap">
    <PeriodicTable
      heatmap_values={mp_elem_counts}
      color_scale={color_scale[0]}
      {log}
      bind:active_element={active_mp_elem}
    >
      <TableInset slot="inset">
        <label for="log">Log color scale<Toggle id="log" bind:checked={log} /></label>
        <PtableInset element={active_mp_elem} elem_counts={mp_elem_counts} />
        <ColorBar
          text_side="top"
          color_scale={color_scale[0]}
          tick_labels={5}
          precision={3}
          range={[0, Math.max(...Object.values(mp_elem_counts))]}
          style="width: 85%; margin: 0 2em 2em;"
        />
      </TableInset>
    </PeriodicTable>
    <p>
      The difference between count modes is best explained by example. <code
        >occurrence</code
      >
      mode maps Fe2O3 to (Fe: 1, O: 1), <code>composition</code> mode maps it to (Fe: 2, O:
      3).
    </p>
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

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
    justify-items: center;
  }
</style>
