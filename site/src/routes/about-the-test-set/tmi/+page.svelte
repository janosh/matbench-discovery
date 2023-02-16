<script lang="ts">
  import { ElemCountInset } from '$lib'
  import {
    ColorScaleSelect,
    PeriodicTable,
    TableInset,
    type ChemicalElement,
  } from 'elementari'
  import { RadioButtons, Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'

  const elem_counts: Record<string, number[]> = {}
  for (const [path, json] of Object.entries(
    import.meta.glob(`../wbm-element-counts-*=*.json`, { eager: true, as: `raw` })
  )) {
    const split = path.split(`.json`)[0].split(`-`).at(-1) as string
    if (!split || !split?.includes(`=`)) console.error(`Invalid path: ${path}`)
    elem_counts[split] = Object.values(JSON.parse(json))
  }

  let arity_keys = Object.keys(elem_counts).filter((k) => k.startsWith(`arity=`))
  let batch_keys = Object.keys(elem_counts).filter((k) => k.startsWith(`batch=`))
  let log = false // log color scale
  let filter = arity_keys[0]
  let selected = [`Turbo`]
  let active_element: ChemicalElement
  $: color_scale = selected[0]
  $: active_counts = elem_counts[filter]

  const style = `display: flex; gap: 5pt; place-items: center; place-content: center;`

  export const snapshot: Snapshot = {
    capture: () => ({ filter, log }),
    restore: (values) => ({ log, filter } = values),
  }
  const radio_style = `display: inline-flex; margin: 2pt 8pt; border-radius: 3pt;`
</script>

<h1>Too Much Information</h1>

Stuff that didn't make the cut into the main page describing the WBM test set.

<h2>WBM Element Counts for <code>{filter}</code></h2>

Filter WBM element counts by composition arity (how many elements in the formula) or batch
index (which iteration of elemental substitution the structure was generated in).

<ColorScaleSelect bind:selected />
<ul>
  <li>
    composition arity:
    <RadioButtons style={radio_style} options={arity_keys} bind:selected={filter}>
      <strong slot="option" let:value let:active class:active>
        {value.split(`=`)[1]}</strong
      >
    </RadioButtons>
  </li>
  <li>
    batch index
    <RadioButtons style={radio_style} options={batch_keys} bind:selected={filter}>
      <strong slot="option" let:value let:active class:active>
        {value.split(`=`)[1]}</strong
      >
    </RadioButtons>
  </li>
</ul>

<PeriodicTable heatmap_values={active_counts} {log} {color_scale} bind:active_element>
  <TableInset slot="inset">
    <ElemCountInset element={active_element} elem_counts={active_counts} />
    <span {style}>Log color scale<Toggle bind:checked={log} /></span>
  </TableInset>
</PeriodicTable>

<style>
  span {
    display: flex;
    gap: 1ex;
  }
  strong {
    background-color: rgba(255, 255, 255, 0.1);
    padding: 3pt 4pt;
  }
  strong.active {
    background-color: teal;
  }
</style>
