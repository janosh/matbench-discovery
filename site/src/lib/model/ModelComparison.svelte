<script lang="ts">
  import { page } from '$app/state'
  import { ALL_METRICS } from '$lib/labels'
  import { is_finite_num } from '$lib/metrics'
  import {
    COMPARE_GROUPS,
    COST_ROWS,
    compare_cells,
    comparison,
    row_value,
    type CompareRow,
  } from '$lib/model-comparison.svelte'
  import ModelSelect from '$lib/ModelSelect.svelte'
  import { ACTIVE_MODELS, MODELS } from '$lib/models.svelte'
  import { rank_color, RANKED_METRICS } from '$lib/rankings'
  import { pareto_staircase } from '$lib/sota'
  import type { ModelData } from '$lib/types'
  import { ScatterPlot, strip_html } from 'matterviz'
  import type { DataSeries, InternalPoint } from 'matterviz/plot'
  import { Icon, Sheet } from 'svelte-widgets'
  import { tooltip } from 'svelte-widgets/attachments'
  import { Cross, Scale } from 'svelte-widgets/icons'

  let models = $derived(comparison.models)
  let n_models = $derived(models.length)
  const has_value = (row: CompareRow, model: ModelData) => {
    const value = row_value(row, model)
    return is_finite_num(value) || (typeof value === `string` && value !== ``)
  }
  // hide rows no compared model has a value for, then groups left empty
  let groups = $derived(
    COMPARE_GROUPS.map((group) => ({
      ...group,
      rows: group.rows.filter((row) => models.some((model) => has_value(row, model))),
    })).filter((group) => group.rows.length > 0),
  )

  // MultiSelect matches options by value, so chips may be rebuilt from the selection
  const option_for = ({ model_key, model_name, lifecycle }: ModelData) => ({
    value: model_key,
    label: `${model_name}${lifecycle === `active` ? `` : ` (${lifecycle})`}`,
  })
  const model_options = MODELS.map(option_for)
  let selected_options = $derived(models.map(option_for))

  // --- cost vs accuracy scatter: whole leaderboard in grey, compared models on top ---
  const all_rows = COMPARE_GROUPS.flatMap((group) => group.rows)
  const row_by_key = Object.fromEntries(all_rows.map((row) => [row.key, row]))
  const accuracy_rows = all_rows.filter((row) => row.better && !COST_ROWS.includes(row))
  const axis_label = ({ label, unit }: CompareRow) =>
    `${label}${unit ? ` (${unit})` : ``}`
  let x_key = $state(COST_ROWS[0].key)
  let y_key = $state<string>()
  // until the user picks one, the y-axis follows the task page the drawer was opened on
  let page_metric = $derived(
    RANKED_METRICS.find((metric) => metric.rank_href === page.url.pathname)?.key ??
      ALL_METRICS.CPS.key,
  )
  let x_row = $derived(row_by_key[x_key])
  let y_row = $derived(row_by_key[y_key ?? page_metric])

  // compared models last so they paint over the field
  let points = $derived(
    [...new Set([...ACTIVE_MODELS, ...models])]
      .flatMap((model) => {
        const [x, y] = [row_value(x_row, model), row_value(y_row, model)]
        if (!is_finite_num(x) || !is_finite_num(y)) return []
        return [{ x, y, model, compared: comparison.keys.has(model.model_key) }]
      })
      .toSorted((pt1, pt2) => Number(pt1.compared) - Number(pt2.compared)),
  )
  // log scale when all positive and spanning at least two decades (params, wall times)
  const log_scale = (values: number[]) => {
    const [min, max] = [Math.min(...values), Math.max(...values)]
    return min > 0 && max >= 100 * min ? `log` : `linear`
  }
  let series: DataSeries<ModelData>[] = $derived.by(() => {
    const field: DataSeries<ModelData> = {
      x: points.map((pt) => pt.x),
      y: points.map((pt) => pt.y),
      markers: `points`,
      metadata: points.map((pt) => pt.model),
      point_style: points.map(({ model, compared }) =>
        compared
          ? { fill: model.color, radius: 7, stroke: `white`, stroke_width: 1.5 }
          : { fill: `gray`, radius: 3.5, fill_opacity: 0.4 },
      ),
      point_label: points.map(({ model, compared }) =>
        compared
          ? { text: model.model_name, font_size: `12px`, auto_placement: true }
          : {},
      ),
    }
    const frontier =
      x_row.better && y_row.better && pareto_staircase(points, x_row.better, y_row.better)
    if (!frontier) return [field]
    const line_style = { stroke: `gray`, stroke_width: 1, line_dash: `4 3` }
    return [field, { ...frontier, markers: `line`, line_style }]
  })
  // frontier corners aren't models: don't snap the tooltip to them
  let tooltip_point: InternalPoint<ModelData> | null = $state(null)
  $effect(() => {
    if (tooltip_point && !tooltip_point.metadata) tooltip_point = null
  })
</script>

<Sheet
  bind:open={comparison.open}
  aria-label="Model comparison"
  style="--sheet-size: min(68rem, 100vw); --sheet-bg: var(--page-bg); --sheet-section-padding: 0.75rem 1rem; --sheet-content-padding: 0 1rem 1rem"
>
  {#snippet trigger(trigger_props)}
    {#if n_models > 0 && !comparison.open}
      <div class="launcher">
        <button {...trigger_props}>
          <Icon icon={Scale} /> Compare {n_models} model{n_models === 1 ? `` : `s`}
        </button>
        <button
          aria-label="Clear model comparison"
          onclick={() => comparison.keys.clear()}
        >
          <Icon icon={Cross} />
        </button>
      </div>
    {/if}
  {/snippet}
  {#snippet header({ close })}
    <div class="head">
      <h2>Compare models</h2>
      <ModelSelect
        options={model_options}
        placeholder="Add models…"
        style="flex: 1; min-width: min(20rem, 100%); border: 1px solid var(--border)"
        bind:selected={
          () => selected_options,
          (options: typeof model_options) =>
            comparison.set(options.map((opt) => opt.value))
        }
      />
      <button aria-label="Close model comparison" onclick={close}>
        <Icon icon={Cross} />
      </button>
    </div>
  {/snippet}

  {#if n_models === 0}
    <p class="empty">
      Pick models above, double-click rows in any leaderboard table, or use the Compare
      button on model pages.
    </p>
  {:else}
    <table>
      <thead>
        <tr>
          <th></th>
          {#each models as model (model.model_key)}
            <th>
              <span style:color={model.color} aria-hidden="true">●</span>
              <a href="/models/{model.model_key}">{model.model_name}</a>
              <button
                aria-label="Remove {model.model_name} from comparison"
                onclick={() => comparison.toggle(model.model_key)}
              >
                <Icon icon={Cross} />
              </button>
            </th>
          {/each}
        </tr>
      </thead>
      {#each groups as group (group.title)}
        <tbody>
          <tr class="group">
            <th colspan={n_models + 1}>
              {#if group.href}<a href={group.href}>{group.title}</a
                >{:else}{group.title}{/if}
            </th>
          </tr>
          {#each group.rows as row (row.key)}
            <tr>
              <th title={row.description} {@attach tooltip({ allow_html: true })}>
                {@html row.label}{#if row.unit}<small> {row.unit}</small>{/if}
              </th>
              {#each compare_cells(row, models) as cell, idx (models[idx].model_key)}
                <td class:best={cell.best} class:text={!row.better}>
                  {cell.text}
                  {#if cell.rank && cell.n}
                    <small
                      style:color={rank_color(cell.rank, cell.n)}
                      title="Ranked {cell.rank} of {cell.n} leaderboard models"
                    >
                      #{cell.rank}
                    </small>
                  {/if}
                </td>
              {/each}
            </tr>
          {/each}
        </tbody>
      {/each}
    </table>

    <h3>Cost vs accuracy</h3>
    <p>
      Every leaderboard model in grey, compared models highlighted; the dashed line traces
      the Pareto frontier. Click a point to add or remove its model.
    </p>
    <!-- plain selects: matterviz's in-plot axis menus portal to <body>, which a modal
    dialog renders inert -->
    <label>
      X axis
      <select value={x_key} onchange={(event) => (x_key = event.currentTarget.value)}>
        {#each COST_ROWS as row (row.key)}
          <option value={row.key}>{strip_html(axis_label(row))}</option>
        {/each}
      </select>
    </label>
    <label>
      Y axis
      <select value={y_row.key} onchange={(event) => (y_key = event.currentTarget.value)}>
        {#each accuracy_rows as row (row.key)}
          <option value={row.key}>{strip_html(axis_label(row))}</option>
        {/each}
      </select>
    </label>
    <ScatterPlot
      {series}
      bind:tooltip_point
      style="height: 360px"
      show_controls={false}
      legend={null}
      x_axis={{
        label: axis_label(x_row),
        format: x_row.format,
        scale_type: log_scale(points.map((pt) => pt.x)),
      }}
      y_axis={{ label: axis_label(y_row), format: y_row.format }}
      point_events={{
        onclick: ({ point }) => {
          const key = point.metadata?.model_key
          if (typeof key === `string`) comparison.toggle(key)
        },
      }}
    >
      {#snippet tooltip({ x_formatted, y_formatted, metadata })}
        {#if metadata}
          <strong>{metadata.model_name}</strong><br />
          {@html x_row.label}: {x_formatted}<br />
          {@html y_row.label}: {y_formatted}<br />
          {@const verb = comparison.keys.has(String(metadata.model_key))
            ? `remove`
            : `add`}
          <small>click to {verb}</small>
        {/if}
      {/snippet}
    </ScatterPlot>
  {/if}
</Sheet>

<style>
  .launcher {
    position: fixed;
    z-index: 20;
    inset: auto 1rem 1rem auto;
    display: flex;
    align-items: center;
    border: 1px solid var(--link-color);
    border-radius: 2em;
    background: var(--page-bg);
    box-shadow: 0 4px 18px var(--shadow);
  }
  .launcher button {
    display: inline-flex;
    align-items: center;
    gap: 0.4em;
    padding: 0.55em 0.9em;
    color: var(--link-color);
    font-weight: 600;
  }
  .launcher button + button {
    padding-left: 0;
  }
  .head {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.5rem 1rem;
  }
  .head h2 {
    margin: 0;
    font-size: 1.2rem;
  }
  /* icon-only buttons (close, clear, remove model) sit flat on the page background */
  button {
    background: transparent;
  }
  :is(.head, thead) button {
    padding: 3pt;
    vertical-align: middle;
  }
  thead button {
    opacity: 0.6;
  }
  .empty {
    margin: 3rem auto;
    max-width: 30rem;
    text-align: center;
    color: var(--text-secondary);
  }
  /* undo the global block+scroll table reset: the sheet body is the scroll container so
     the sticky header row and metric column stay put while scrolling */
  table {
    display: table;
    overflow: visible;
    width: 100%;
    margin: 1rem 0 2rem;
    font-size: 0.9em;
  }
  tbody tr {
    background: none; /* no global zebra stripes: sticky cells need an opaque bg */
  }
  th,
  td {
    padding: 3pt 8pt;
    border: 0;
    border-bottom: 1px solid var(--border);
    text-align: left;
    vertical-align: top;
    white-space: nowrap;
    font-variant-numeric: tabular-nums;
  }
  th:first-child {
    position: sticky;
    left: 0;
    background: var(--page-bg);
  }
  thead th {
    position: sticky;
    top: 0;
    z-index: 1;
    background: var(--page-bg);
    font-size: 1.05em;
  }
  thead th:first-child {
    z-index: 2;
  }
  tbody th {
    font-weight: normal;
    color: var(--text-secondary);
  }
  tr.group th {
    padding-top: 1.2em;
    font-weight: 600;
    color: var(--text-color);
  }
  td.text {
    white-space: normal;
    min-width: 11rem;
  }
  td.best {
    font-weight: 600;
    background: color-mix(in srgb, var(--link-color) 12%, transparent);
  }
  td small {
    margin-left: 4pt;
    font-size: 0.75em;
  }
  h3 {
    margin: 0 0 0.3em;
  }
  h3 + p {
    margin: 0 0 0.5em;
    font-size: 0.9em;
    color: var(--text-secondary);
  }
  label {
    margin-right: 1.5em;
    font-size: 0.9em;
  }
  select {
    margin-left: 0.3em;
  }
</style>
