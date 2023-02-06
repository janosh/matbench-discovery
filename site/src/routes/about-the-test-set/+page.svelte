<script lang="ts">
  import { browser } from '$app/environment'
  import FormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import DataReadme from '$root/data/wbm/readme.md'
  import WbmEachHist from '$site/src/figs/wbm-each-hist.svelte'
  import { Toggle } from 'svelte-zoo'
  import type { ChemicalElement } from 'sveriodic-table'
  import { ColorScaleSelect, PeriodicTable, TableInset } from 'sveriodic-table'
  import { pretty_num } from 'sveriodic-table/labels'
  import mp_elem_counts from './2023-01-08-mp-element-counts.json'
  import wbm_elem_counts from './2023-01-08-wbm-element-counts.json'

  let log = true // log color scale
  let active_mp_elem: ChemicalElement
  let active_wbm_elem: ChemicalElement

  const wbm_heat_vals: number[] = Object.values(wbm_elem_counts)
  const mp_heat_vals: number[] = Object.values(mp_elem_counts)
  let color_scale: string
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if browser}
      <FormEnergyHist />
    {/if}
  </svelte:fragment>
  <svelte:fragment slot="wbm-elements-heatmap">
    <span>Log color scale <Toggle bind:checked={log} /></span>
    <PeriodicTable
      heatmap_values={wbm_heat_vals}
      {color_scale}
      {log}
      bind:active_element={active_wbm_elem}
    >
      <TableInset slot="inset" grid_row="3">
        {#if active_wbm_elem?.name}
          <strong>
            {active_wbm_elem?.name}: {pretty_num(
              wbm_elem_counts[active_wbm_elem?.symbol]
            )}
            <!-- compute percent of total -->
            {#if wbm_elem_counts[active_wbm_elem?.symbol] > 0}
              {@const total = wbm_heat_vals.reduce((a, b) => a + b, 0)}
              ({pretty_num((wbm_elem_counts[active_wbm_elem?.symbol] / total) * 100)}%)
            {/if}
          </strong>
        {/if}
      </TableInset>
    </PeriodicTable>
    <ColorScaleSelect bind:value={color_scale} selected={[`Turbo`]} />
  </svelte:fragment>
  <svelte:fragment slot="mp-elements-heatmap">
    <span>Log color scale <Toggle bind:checked={log} /></span>
    <PeriodicTable
      heatmap_values={mp_heat_vals}
      {color_scale}
      {log}
      bind:active_element={active_mp_elem}
    >
      <TableInset slot="inset" grid_row="3">
        {#if active_mp_elem?.name}
          <strong>
            {active_mp_elem?.name}: {pretty_num(mp_elem_counts[active_mp_elem?.symbol])}
            <!-- compute percent of total -->
            {#if mp_elem_counts[active_mp_elem?.symbol] > 0}
              {@const total = wbm_heat_vals.reduce((a, b) => a + b, 0)}
              ({pretty_num((mp_elem_counts[active_mp_elem?.symbol] / total) * 100)}%)
            {/if}
          </strong>
        {/if}
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
  span {
    display: flex;
    gap: 1ex;
    position: absolute;
    left: 50%;
    transform: translateX(-50%) translateY(100%);
    z-index: 1;
  }
  strong {
    text-align: center;
  }
</style>
