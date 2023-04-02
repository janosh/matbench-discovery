<script lang="ts">
  import { browser } from '$app/environment'
  import BarElementCounts from '$figs/bar-element-counts-mp+wbm-normalized=False.svelte'
  import BarElementCountsNormalized from '$figs/bar-element-counts-mp+wbm-normalized=True.svelte'
  import { PtableInset } from '$lib'
  import type { ChemicalElement } from 'elementari'
  import { ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import { RadioButtons, Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'

  const elem_counts: Record<string, number[]> = {}
  for (const [path, json] of Object.entries(
    import.meta.glob(`../wbm-element-counts-*=*.json`, { eager: true, as: `raw` })
  )) {
    const split = path?.split(`.json`)[0]?.split(`-`).at(-1) as string
    if (!split || !split?.includes(`=`)) console.error(`Invalid path: ${path}`)
    elem_counts[split] = JSON.parse(json)
  }

  let arity_keys = Object.keys(elem_counts).filter((k) => k.startsWith(`arity=`))
  let batch_keys = Object.keys(elem_counts).filter((k) => k.startsWith(`batch=`))
  let log = false // log color scale
  let filter = arity_keys[0]
  let color_scale = [`Inferno`]
  let active_element: ChemicalElement
  $: active_counts = elem_counts[filter]
  let normalized_bar_counts: boolean = false

  const style = `display: flex; place-items: center; place-content: center;`

  export const snapshot: Snapshot = {
    capture: () => ({ filter, log, color_scale }),
    restore: (values) => ({ filter, log, color_scale } = values),
  }
  const radio_style = `display: inline-flex; margin: 2pt 8pt; border-radius: 3pt;`
</script>

<h1>Too Much Information</h1>

Stuff that didn't make the cut into the main page describing the WBM test set.

<h2>WBM Element Counts for <code>{filter}</code></h2>

Filter WBM element counts by composition<strong>arity</strong> (how many elements in the
formula) or <strong>batch index</strong> (which iteration of elemental substitution the
structure was generated in).

<ColorScaleSelect bind:selected={color_scale} />
<ul>
  <li>
    composition arity:
    <RadioButtons style={radio_style} options={arity_keys} bind:selected={filter}>
      <strong slot="option" let:option let:active class:active>
        {option?.split(`=`)[1]}</strong
      >
    </RadioButtons>
  </li>
  <li>
    batch index
    <RadioButtons style={radio_style} options={batch_keys} bind:selected={filter}>
      <strong slot="option" let:option let:active class:active>
        {option?.split(`=`)[1]}</strong
      >
    </RadioButtons>
  </li>
</ul>

<PeriodicTable
  heatmap_values={active_counts}
  {log}
  color_scale={color_scale[0]}
  bind:active_element
>
  <TableInset slot="inset">
    <PtableInset element={active_element} elem_counts={active_counts} />
    <span {style}>Log color scale<Toggle bind:checked={log} /></span>
  </TableInset>
</PeriodicTable>

<h2>Element Counts</h2>

<label>
  Normalize by data set size
  <input type="checkbox" bind:checked={normalized_bar_counts} />
</label>

{#if browser}
  {#if normalized_bar_counts}
    <BarElementCountsNormalized />
  {:else}
    <BarElementCounts />
  {/if}
{/if}

<style>
  span {
    display: flex;
    gap: 1ex;
  }
  ul > li strong {
    background-color: rgba(255, 255, 255, 0.1);
    padding: 3pt 4pt;
  }
  ul > li strong.active {
    background-color: teal;
  }
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    place-items: center;
  }
</style>
