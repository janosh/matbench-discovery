<script lang="ts">
  import elem_counts_bar from '$figs/element-counts-mp-vs-wbm.json.gz'
  import { PtableInset } from '$lib'
  import type { ChemicalElement } from 'matterviz'
  import { ColorScaleSelect, PeriodicTable, TableInset } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { BarPlot } from 'matterviz/plot'
  import { Toggle } from 'svelte-multiselect'
  import {
    bind_url_params,
    bool_from_param,
    bool_url_entry,
    url_color_scale,
    valid_query_param,
  } from '$lib/url-state.svelte'

  const raw_counts = import.meta.glob<Record<string, number>>(
    `../wbm-element-counts-*=*.json`,
    { eager: true, import: 'default' },
  )
  // key files by their `arity=N`/`batch=N` filename suffix
  const elem_counts: Record<string, Record<string, number>> = $state({})
  for (const [key, value] of Object.entries(raw_counts)) {
    const short_key = key.split(`-`).at(-1)?.split(`.`)[0]
    if (short_key) elem_counts[short_key] = value
  }

  const all_keys = Object.keys(elem_counts)
  const starts_with_arity = (key: string) => key.startsWith(`arity=`)
  let arity_keys = all_keys.filter(starts_with_arity)
  let batch_keys = all_keys.filter((key) => key.startsWith(`batch=`))
  const default_filter = all_keys.find(starts_with_arity) ?? all_keys[0] ?? ``
  let log = $state(false) // Log color scale
  let filter = $state(default_filter)
  let color_scale = $state<D3InterpolateName>(url_color_scale.default)
  let active_element: ChemicalElement | null = $state(null)
  let active_counts = $derived(elem_counts[filter])
  let normalized_bar_counts: boolean = $state(false)

  const read_url_params = (params: URLSearchParams) => {
    filter = valid_query_param(params, `filter`, default_filter, new Set(all_keys))
    log = bool_from_param(params, `log`)
    normalized_bar_counts = bool_from_param(params, `normalized`)
    color_scale = url_color_scale.read(params)
  }
  bind_url_params(read_url_params, () => [
    [`filter`, filter, default_filter],
    bool_url_entry(`log`, log),
    bool_url_entry(`normalized`, normalized_bar_counts),
    url_color_scale.entry(color_scale),
  ])

  const style = `display: flex; place-items: center; place-content: center;`
</script>

<h1>Too Much Information</h1>

<p>Stuff that didn't make the cut into the main page describing the WBM test set.</p>

<h2>WBM Element Counts for <code>{filter}</code></h2>

<p>
  Filter WBM element counts by composition<strong>arity</strong> (how many elements in the
  formula) or <strong>batch index</strong> (which iteration of elemental substitution the structure
  was generated in).
</p>

<ColorScaleSelect bind:value={color_scale} selected={[color_scale]} />

<form>
  <span>
    number of elements
    {#each arity_keys as value (value)}
      <label>
        <input type="radio" name="filter" {value} bind:group={filter} />
        <strong class:active={filter === value}>{value.split(`=`).at(-1)}</strong>
      </label>
    {/each}
  </span>
  <span>
    batch index
    {#each batch_keys as value (value)}
      <label>
        <input type="radio" name="filter" {value} bind:group={filter} />
        <strong class:active={filter === value}>{value.split(`=`).at(-1)}</strong>
      </label>
    {/each}
  </span>
</form>

<PeriodicTable
  heatmap_values={active_counts}
  {log}
  {color_scale}
  bind:active_element
  show_photo={false}
  missing_color="rgba(255,255,255,0.3)"
>
  {#snippet inset()}
    <TableInset>
      {#if active_element}
        <PtableInset element={active_element} elem_counts={active_counts} />
      {/if}
      <span {style}>Log color scale<Toggle bind:checked={log} /></span>
    </TableInset>
  {/snippet}
</PeriodicTable>

<h2>Element Counts</h2>

<label>
  Normalize by data set size
  <input type="checkbox" bind:checked={normalized_bar_counts} />
</label>

<BarPlot
  series={elem_counts_bar[normalized_bar_counts ? `normalized` : `raw`]}
  mode="grouped"
  x_axis={{ label: `Element`, range: [0, null] }}
  y_axis={{
    label: normalized_bar_counts ? `Share of Structures (%)` : `Number of Structures`,
    format: `~s`,
  }}
  show_legend
  show_controls={false}
  style="height: 360px"
/>

<style>
  span {
    display: flex;
    gap: 1ex;
    place-items: center;
  }
  form {
    display: flex;
    gap: 5cqw;
    place-content: center;
  }
  form > span strong.active {
    background-color: var(--btn-bg);
  }
  label {
    display: flex;
    gap: 4pt;
  }
</style>
