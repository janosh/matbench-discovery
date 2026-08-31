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
  import { ACTIVE_MODELS, MODELS } from '$lib/models.svelte'
  import { rank_color, RANKED_METRICS } from '$lib/rankings'
  import { pareto_staircase } from '$lib/sota'
  import type { ModelData } from '$lib/types'
  import { format_num, ScatterPlot, strip_html } from 'matterviz'
  import type { DataSeries, InternalPoint } from 'matterviz/plot'
  import { tick, untrack } from 'svelte'
  import { Dialog, Icon, MultiSelect, type Option } from 'svelte-widgets'
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
  const option_for = ({ model_key, model_name, lifecycle, CPS, dates }: ModelData) => ({
    value: model_key,
    label: `${model_name}${lifecycle === `active` ? `` : ` (${lifecycle})`}`,
    cps: CPS,
    added: dates.benchmark_added ?? ``,
  })
  type ModelOption = ReturnType<typeof option_for>
  // dropdown order: best CPS first (unscored models last) or most recently added first
  let option_sort = $state<`cps` | `added`>(`cps`)
  const by_cps = (opt_1: ModelOption, opt_2: ModelOption) =>
    (is_finite_num(opt_2.cps) ? opt_2.cps : -Infinity) -
    (is_finite_num(opt_1.cps) ? opt_1.cps : -Infinity)
  const by_added = (opt_1: ModelOption, opt_2: ModelOption) =>
    opt_2.added.localeCompare(opt_1.added)
  // CPS is live-reweighted, so options are derived rather than built once
  let model_options = $derived(
    MODELS.map(option_for).toSorted(option_sort === `cps` ? by_cps : by_added),
  )
  let selected_options = $derived(models.map(option_for))
  const external = { target: `_blank`, rel: `noopener` }
  // opened with nothing to compare yet (empty or a single model), the picker is the next
  // step, so it takes focus (which also drops its option list open). n_models is untracked
  // so removing models while the dialog is open doesn't yank focus to the picker
  let picker_input = $state<HTMLInputElement | null>(null)
  $effect(() => {
    if (comparison.open && untrack(() => n_models) <= 1) {
      void tick().then(() => picker_input?.focus())
    }
  })
  const MAX_TRAY_CHIPS = 4

  // --- cost vs accuracy scatter: whole leaderboard in grey, compared models on top ---
  const all_rows = COMPARE_GROUPS.flatMap((group) => group.rows)
  const row_by_key = Object.fromEntries(all_rows.map((row) => [row.key, row]))
  const accuracy_rows = all_rows.filter((row) => row.better && !COST_ROWS.includes(row))
  const axis_label = ({ label, unit }: CompareRow) =>
    `${label}${unit ? ` (${unit})` : ``}`
  let x_key = $state(COST_ROWS[0].key)
  let y_key = $state<string>()
  // until the user picks one, the y-axis follows the task page the dialog was opened on
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

<Dialog
  bind:open={comparison.open}
  aria-label="Model comparison"
  style="--dialog-width: min(68rem, calc(100vw - 2rem)); --dialog-radius: 12px; --dialog-bg: var(--page-bg); --dialog-section-padding: 0.5rem 1rem 0; --dialog-content-padding: 0 1rem 1rem"
>
  <!-- selection tray: shows up with the first pick so the next step is never a guess -->
  {#snippet trigger(trigger_props)}
    {#if n_models > 0 && !comparison.open}
      <div class="tray">
        <!-- text-only live region: a status role on the whole tray would re-announce every
        chip and button label on each add/remove -->
        <span class="visually-hidden" role="status">
          {n_models}
          {n_models === 1 ? `model` : `models`} selected for comparison
        </span>
        <ul>
          {#each models.slice(0, MAX_TRAY_CHIPS) as model (model.model_key)}
            <li>
              <span style:color={model.color} aria-hidden="true">●</span>
              {model.model_name}
              <button
                aria-label="Remove {model.model_name} from comparison"
                onclick={() => comparison.toggle(model.model_key)}
              >
                <Icon icon={Cross} />
              </button>
            </li>
          {/each}
          {#if n_models > MAX_TRAY_CHIPS}
            <li>+{n_models - MAX_TRAY_CHIPS} more</li>
          {/if}
        </ul>
        {#if n_models === 1}
          <span class="hint">Pick a 2nd model to compare against</span>
        {/if}
        <button class="open" {...trigger_props}>
          <Icon icon={Scale} />
          {n_models === 1 ? `Open comparison` : `Compare ${n_models} models →`}
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
      <MultiSelect
        options={model_options}
        placeholder="Add models…"
        style="flex: 1; min-width: min(20rem, 100%); border: 1px solid var(--border)"
        bind:input={picker_input}
        bind:selected={
          () => selected_options,
          (options: ModelOption[]) => comparison.set(options.map((opt) => opt.value))
        }
      >
        {#snippet children({ option, type }: { option: Option; type: string })}
          <!-- MultiSelect isn't generic over our option shape, only ever fed model_options -->
          {@const opt = option as ModelOption}
          {opt.label}
          {#if type === `option`}
            <small class="option-detail">
              {#if is_finite_num(opt.cps)}{format_num(opt.cps, `.3f`)}{/if}
              {#if option_sort === `added`}· {opt.added}{/if}
            </small>
          {/if}
        {/snippet}
      </MultiSelect>
      <label class="option-sort">
        Sort
        <select
          value={option_sort}
          onchange={(event) =>
            (option_sort = event.currentTarget.value === `added` ? `added` : `cps`)}
        >
          <option value="cps">by CPS</option>
          <option value="added">by newest</option>
        </select>
      </label>
      <button aria-label="Close model comparison" onclick={close}>
        <Icon icon={Cross} />
      </button>
    </div>
  {/snippet}

  {#if n_models === 0}
    <p class="empty">
      Pick models above. In any leaderboard table you can also hover a row and click its
      <Icon icon={Scale} style="vertical-align: middle" /> button, right-click the row, or double-click
      it.
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
                  {#if cell.parts}
                    {#each cell.parts as part, part_idx (part_idx)}
                      {#if part.href}
                        <a
                          href={part.href}
                          title={part.title}
                          {...part.href.startsWith(`http`) ? external : {}}
                          {@attach tooltip()}>{part.text}</a
                        >
                      {:else if part.title}
                        <span title={part.title} {@attach tooltip()}>{part.text}</span>
                      {:else}{part.text}{/if}
                    {/each}
                  {:else if cell.title}
                    <span title={cell.title} {@attach tooltip({ allow_html: true })}>
                      {cell.text}
                    </span>
                  {:else}
                    {cell.text}
                  {/if}
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
      y_axis={{
        label: axis_label(y_row),
        format: y_row.format,
        scale_type: log_scale(points.map((pt) => pt.y)),
      }}
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
</Dialog>

<style>
  .tray {
    position: fixed;
    z-index: 20;
    inset: auto 0 1rem;
    width: fit-content;
    max-width: calc(100vw - 2rem);
    margin-inline: auto;
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
    gap: 0.3em 0.6em;
    padding: 0.4em 0.5em 0.4em 0.7em;
    border: 1px solid var(--link-color);
    border-radius: 2em;
    background: var(--page-bg);
    box-shadow: 0 4px 18px var(--shadow);
    font-size: 0.95em;
    transition:
      translate 0.25s,
      opacity 0.25s;
    @starting-style {
      translate: 0 1rem;
      opacity: 0;
    }
  }
  .tray ul {
    display: contents;
  }
  .tray li {
    display: inline-flex;
    align-items: center;
    gap: 0.3em;
    padding: 0.1em 0.2em 0.1em 0.5em;
    border-radius: 1em;
    background: var(--chip-bg);
    white-space: nowrap;
  }
  .tray .hint {
    color: var(--text-secondary);
  }
  .tray .visually-hidden {
    position: absolute;
    width: 1px;
    height: 1px;
    overflow: hidden;
    clip-path: inset(50%);
  }
  .tray button.open {
    display: inline-flex;
    align-items: center;
    gap: 0.4em;
    padding: 0.4em 0.9em;
    border-radius: 2em;
    background: var(--link-color);
    color: var(--page-bg);
    font-weight: 600;
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
  :is(.head, thead, .tray) button {
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
  /* undo the global block+scroll table reset: the dialog body is the scroll container so
     the sticky header row and metric column stay put while scrolling */
  table {
    display: table;
    overflow: visible;
    width: 100%;
    margin: 0 0 2rem;
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
  tbody:first-of-type tr.group th {
    padding-top: 0.4em;
  }
  .option-detail {
    margin-left: 0.5em;
    color: var(--text-secondary);
    font-variant-numeric: tabular-nums;
  }
  .option-sort {
    margin: 0;
    font-size: 0.85em;
    color: var(--text-secondary);
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
