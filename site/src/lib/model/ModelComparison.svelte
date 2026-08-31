<script lang="ts">
  import { page } from '$app/state'
  import { ALL_METRICS, HYPERPARAMS, scatter_options } from '$lib/labels'
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
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import { rank_color, RANKED_METRICS } from '$lib/rankings'
  import type { ModelData } from '$lib/types'
  import { format_num } from 'matterviz'
  import { tick, untrack } from 'svelte'
  import { Dialog, Icon, MultiSelect, type Option } from 'svelte-widgets'
  import { tooltip } from 'svelte-widgets/attachments'
  import { Cross } from 'svelte-widgets/icons'

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

  // --- cost vs accuracy scatter: whole leaderboard dimmed, compared models on top ---
  // cost rows lead the axis menus and carry lower=better so the Pareto frontier can be
  // traced against them; the rest of the usual scatter options follow
  const cost_keys = new Set(COST_ROWS.map((row) => row.key))
  const scatter_axis_options = [
    ...COST_ROWS,
    ...scatter_options.filter((opt) => !cost_keys.has(opt.key)),
  ]
  let x_key = $state(HYPERPARAMS.model_params.key)
  let y_key = $state<string>()
  // until the user picks one, the y-axis follows the task page the dialog was opened on
  let page_metric = $derived(
    RANKED_METRICS.find((metric) => metric.rank_href === page.url.pathname)?.key ??
      ALL_METRICS.CPS.key,
  )
</script>

<Dialog
  bind:open={comparison.open}
  aria-label="Model comparison"
  style="--dialog-width: min(68rem, calc(100vw - 2rem)); --dialog-section-padding: 0.5rem 1rem 0; --dialog-content-padding: 0 1rem 1rem;"
>
  {#snippet header({ close })}
    <div class="head">
      <h2>Compare models</h2>
      <!-- the option list stays inside the dialog (overflow: hidden, no portal past a modal's
        top layer), so its max height plus the header must fit the dialog's min-height -->
      <MultiSelect
        options={model_options}
        placeholder="Add models…"
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
      Pick models above. In any leaderboard table you can also right-click or double-click
      a row to add its model.
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
      Every leaderboard model dimmed, compared models highlighted; the dashed line traces
      the Pareto frontier. Click a point to add or remove its model.
    </p>
    <DynamicScatter
      models={[...new Set([...ACTIVE_MODELS, ...models])]}
      options={scatter_axis_options}
      highlight_keys={comparison.keys}
      show_pareto_frontier
      bleed={false}
      legend={null}
      style="height: 420px"
      bind:x_key
      bind:y_key={() => y_key ?? page_metric, (key) => (y_key = key)}
      point_events={{
        onclick: ({ point }) => {
          const key = point.metadata?.model_key
          if (typeof key === `string`) comparison.toggle(key)
        },
      }}
    />
  {/if}
</Dialog>

<style>
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
  /* selected chips: smaller plain text, the shading moves onto a circular remove button */
  .head :global(ul.selected > li) {
    font-size: 0.85em;
  }
  .head :global(ul.selected > li button.remove) {
    background: light-dark(rgba(100, 120, 140, 0.15), rgba(120, 170, 255, 0.2));
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
</style>
