<script lang="ts">
  import { browser } from '$app/environment'
  import FormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import SpacegroupSunburstMp from '$figs/spacegroup-sunburst-mp.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import WbmEachHist from '$figs/wbm-each-hist.svelte'
  import { PtableInset } from '$lib'
  import DataReadme from '$root/data/wbm/readme.md'
  import type { ChemicalElement } from 'elementari'
  import { ColorBar, ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import Select from 'svelte-multiselect'
  import { Toggle, Tooltip } from 'svelte-zoo'
  import type { Snapshot } from './$types'

  const elem_counts = import.meta.glob(`./*-element-counts-{occu,comp}*.json`, {
    eager: true,
    import: 'default',
  })

  let log = false // log color scale
  let color_scale = [`Inferno`]
  let active_mp_elem: ChemicalElement
  let active_wbm_elem: ChemicalElement
  const count_mode_ops = [`occurrence`, `composition`]
  let count_mode = [count_mode_ops[0]]

  $: mp_elem_counts = elem_counts[`./mp-element-counts-${count_mode[0]}.json`]
  $: wbm_elem_counts = elem_counts[`./wbm-element-counts-${count_mode[0]}.json`]

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
    <Tooltip
      text="occurrence=(Fe: 1, O: 1), composition: Fe2O3=(Fe: 2, O: 3)"
      style="display: inline-block; transform: translate(10cqw, 5ex);"
    >
      <label for="count-mode">Count Mode</label>
    </Tooltip>
    <Select
      id="count-mode"
      bind:selected={count_mode}
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
  </svelte:fragment>
  <svelte:fragment slot="wbm-each-hist">
    {#if browser}
      <WbmEachHist />
    {/if}
  </svelte:fragment>
  <div
    style="display: flex; gap: 1em; justify-content: space-around;"
    slot="spacegroup-sunbursts"
  >
    {#if browser}
      <SpacegroupSunburstMp />
      <SpacegroupSunburstWbm />
    {/if}
  </div>
</DataReadme>

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
    justify-items: center;
  }
</style>
