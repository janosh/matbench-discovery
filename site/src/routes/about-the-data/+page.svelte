<script lang="ts">
  import { browser } from '$app/environment'
  import FormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import WbmEachHist from '$figs/wbm-each-hist.svelte'
  import { ElemCountInset } from '$lib'
  import DataReadme from '$root/data/wbm/readme.md'
  import type { ChemicalElement } from 'elementari'
  import { ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import { Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'
  import mp_elem_counts from './mp-element-counts.json'
  import wbm_elem_counts from './wbm-element-counts.json'

  let log = true // log color scale
  let color_scale = [`Turbo`]
  let active_mp_elem: ChemicalElement
  let active_wbm_elem: ChemicalElement

  export const snapshot: Snapshot = {
    capture: () => ({ color_scale, log }),
    restore: (values) => ({ color_scale, log } = values),
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
        <ElemCountInset element={active_wbm_elem} elem_counts={wbm_elem_counts} />
      </TableInset>
    </PeriodicTable>
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
        <ElemCountInset element={active_mp_elem} elem_counts={mp_elem_counts} />
      </TableInset>
    </PeriodicTable>
  </svelte:fragment>
  <svelte:fragment slot="wbm-each-hist">
    {#if browser}
      <WbmEachHist />
    {/if}
  </svelte:fragment>
</DataReadme>

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
  }
</style>
