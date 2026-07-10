<script lang="ts">
  import type { Label } from '$lib'
  import {
    BUILTIN_PRESETS,
    delete_user_preset,
    save_user_preset,
    user_presets,
  } from '$lib/filter-presets.svelte'
  import { openness_tooltips } from '$lib/metrics'
  import { make_table_filters, MODELS } from '$lib/models.svelte'
  import {
    DEFAULT_TARGETS_PARAM,
    FS_MODES,
    OPENNESS_OPTIONS,
    parse_targets,
    TARGET_OUTPUTS,
    type TargetOutput,
    type UrlTableFilters,
  } from '$lib/url-state.svelte'
  import { Icon } from 'matterviz'
  import { type Label as MvLabel, ToggleMenu } from 'matterviz/table'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    columns = $bindable([]),
    filters = make_table_filters(),
    show_selected_only = $bindable(false),
    selected_count = 0,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    columns?: Label[]
    filters?: UrlTableFilters
    show_selected_only?: boolean
    selected_count?: number
  } = $props()

  const close_on_outside_click = click_outside({
    callback: (node) => ((node as HTMLDetailsElement).open = false),
  })

  // summaries open their pane directly below, so tooltips go above the button
  const summary_tooltip = tooltip({ placement: `top` })

  const train_modes = [`require`, `exclude`] as const
  // static per-category model tallies shown in the filter panels. Semantics mirror
  // UrlTableFilters.matches: missing openness counts as OSOD, targets are parsed with
  // the same parse_targets, and fs_mode `any` counts every model.
  const openness_counts: Record<string, number> = {}
  const training_counts: Record<string, number> = {}
  const target_counts: Record<string, number> = {}
  const fs_mode_counts: Record<string, number> = { any: MODELS.length }
  for (const model of MODELS) {
    const openness = model.openness ?? `OSOD`
    openness_counts[openness] = (openness_counts[openness] ?? 0) + 1
    for (const dataset of model.training_set as string[]) {
      training_counts[dataset] = (training_counts[dataset] ?? 0) + 1
    }
    const { outputs, fs_mode } = parse_targets(model.targets)
    for (const output of outputs) {
      target_counts[output] = (target_counts[output] ?? 0) + 1
    }
    if (fs_mode) fs_mode_counts[fs_mode] = (fs_mode_counts[fs_mode] ?? 0) + 1
  }
  const n_train = $derived(Object.keys(filters.training).length)
  const n_openness = $derived(filters.openness.length)
  // badge shows the active constraints when they differ from the default (require F)
  const targets_badge = $derived(
    filters.targets_param === DEFAULT_TARGETS_PARAM
      ? ``
      : ` (${filters.targets_param || `all`})`,
  )

  let new_preset_name = $state(``)
  function save_current_filters(event: SubmitEvent) {
    event.preventDefault()
    const name = new_preset_name.trim()
    if (!name) return
    save_user_preset(name, filters.as_preset)
    new_preset_name = ``
  }
</script>

<div class="table-controls" {...rest}>
  {#if selected_count > 0 || show_selected_only}
    <label>
      <input
        type="checkbox"
        bind:checked={show_selected_only}
        aria-label="Toggle between showing only selected models and all models"
      />
      {show_selected_only ? `Show all` : `Show only ${selected_count} selected`}
    </label>
  {/if}

  <details class="filter-menu" {@attach close_on_outside_click}>
    <summary
      title="Filter models by the datasets they were trained on"
      {@attach summary_tooltip}
    >
      Training data{n_train ? ` (${n_train})` : ``}
    </summary>
    <div class="dropdown train-grid">
      <span class="hint">
        <em>require</em> = model's training set must include this dataset,
        <em>exclude</em> = hide models trained on it
      </span>
      {#each train_modes as mode, mode_idx (mode)}
        <span class="col-head" style:grid-column={mode_idx + 2}>{mode}</span>
      {/each}
      {#each filters.training_sets as dataset_key (dataset_key)}
        <span>{dataset_key} ({training_counts[dataset_key] ?? 0})</span>
        {#each train_modes as mode (mode)}
          <input
            type="checkbox"
            aria-label="{mode} {dataset_key}"
            checked={filters.training[dataset_key] === mode}
            onchange={() => filters.set_training(dataset_key, mode)}
          />
        {/each}
      {/each}
    </div>
  </details>

  <details class="filter-menu" {@attach close_on_outside_click}>
    <summary
      title="Filter models by whether their source code and training data are open"
      {@attach summary_tooltip}
    >
      Openness{n_openness < OPENNESS_OPTIONS.length
        ? ` (${n_openness}/${OPENNESS_OPTIONS.length})`
        : ``}
    </summary>
    <div class="dropdown">
      {#each OPENNESS_OPTIONS as openness (openness)}
        <label class="filter-row" title={openness_tooltips[openness]} {@attach tooltip()}>
          <input
            type="checkbox"
            checked={filters.openness.includes(openness)}
            onchange={() => filters.toggle_openness(openness)}
          />
          {openness} ({openness_counts[openness] ?? 0})
        </label>
      {/each}
    </div>
  </details>

  <details class="filter-menu" {@attach close_on_outside_click}>
    <summary
      title="Filter models by which quantities they predict and how forces/stress are computed"
      {@attach summary_tooltip}
    >
      Targets{targets_badge}
    </summary>
    <div class="dropdown train-grid">
      <span class="hint">
        Every model predicts energy (E). <em>require</em>/<em>exclude</em> filter by the other
        predicted outputs; forces are required by default (hides energy-only models)
      </span>
      {#each train_modes as mode, mode_idx (mode)}
        <span class="col-head" style:grid-column={mode_idx + 2}>{mode}</span>
      {/each}
      {#each Object.entries(TARGET_OUTPUTS) as [key, label] (key)}
        <span>{label} ({key}) ({target_counts[key] ?? 0})</span>
        {#each train_modes as mode (mode)}
          <input
            type="checkbox"
            aria-label="{mode} {label}"
            checked={filters.targets[key as TargetOutput] === mode}
            onchange={() => filters.set_target(key as TargetOutput, mode)}
          />
        {/each}
      {/each}
      <span class="fs-mode">
        forces/stress via
        <span class="fs-mode-options">
          {#each FS_MODES as mode (mode)}
            <label>
              {mode} ({fs_mode_counts[mode] ?? 0})
              <input
                type="radio"
                checked={filters.fs_mode === mode}
                onchange={() => (filters.fs_mode = mode)}
              />
            </label>
          {/each}
        </span>
      </span>
    </div>
  </details>

  <details class="filter-menu" {@attach close_on_outside_click}>
    <summary
      title="Apply a saved filter combination or save the current one"
      {@attach summary_tooltip}
    >
      Presets
    </summary>
    <div class="dropdown">
      {#each Object.entries( { ...BUILTIN_PRESETS, ...user_presets }, ) as [name, preset] (name)}
        <span class="filter-row">
          <button
            class="preset"
            onclick={() => filters.apply(preset)}
            title={preset.description}
            {@attach tooltip()}
          >
            {name}
          </button>
          {#if name in user_presets}
            <button
              class="delete-preset"
              aria-label="Delete preset {name}"
              onclick={() => delete_user_preset(name)}
            >
              <Icon icon="Cross" />
            </button>
          {/if}
        </span>
      {/each}
      <form onsubmit={save_current_filters}>
        <input
          placeholder="Save current filters as…"
          aria-label="New preset name"
          bind:value={new_preset_name}
        />
        <button disabled={!new_preset_name.trim()}>Save</button>
      </form>
    </div>
  </details>

  {#if filters.n_active > 0}
    <button
      class="clear-filters"
      onclick={() => filters.clear()}
      title="Reset training-data, openness and target filters"
      {@attach tooltip()}
    >
      <Icon icon="Cross" /> clear filters
    </button>
  {/if}

  <label>
    <input
      type="checkbox"
      bind:checked={filters.show_heatmap}
      aria-label="Toggle heatmap colors"
    />
    Heatmap
  </label>

  <ToggleMenu bind:columns={columns as MvLabel[]} />
</div>

<style>
  .table-controls {
    display: inline-flex;
    flex-wrap: wrap;
    justify-content: end;
    gap: 4pt 12pt;
    align-items: center;
    font-size: clamp(9pt, 1.4cqw, 11pt);
  }
  details.filter-menu {
    position: relative;
  }
  details.filter-menu summary {
    cursor: pointer;
    list-style: none;
    padding: 1pt 6pt;
    border-radius: 4px;
    background: var(--btn-bg);
  }
  details.filter-menu[open] summary {
    background: color-mix(in srgb, var(--link-color) 25%, transparent);
  }
  details.filter-menu .dropdown {
    position: absolute;
    right: 0;
    z-index: 6;
    display: grid;
    gap: 2pt;
    min-width: max-content;
    margin-top: 4px;
    padding: 4pt 6pt;
    background: var(--page-bg);
    border: 1px solid var(--border);
    border-radius: 5px;
    box-shadow: 0 0 10px var(--shadow);
  }
  details.filter-menu .hint {
    font-size: 0.85em;
    opacity: 0.75;
    text-wrap: balance;
    margin-bottom: 1pt;
  }
  /* training-data panel: 3-column grid (dataset | require | exclude) */
  .train-grid {
    grid-template-columns: auto repeat(2, min-content);
    gap: 2pt 8pt;
    align-items: center;
  }
  .train-grid .hint {
    grid-column: 1 / -1;
    /* contribute zero to grid column sizing (the dataset rows set the panel width),
    then stretch to whatever width they produce and wrap */
    width: 0;
    min-width: 100%;
    text-wrap: wrap;
  }
  .train-grid :is(.col-head, input) {
    justify-self: center;
  }
  .train-grid .col-head {
    font-style: italic;
    opacity: 0.75;
  }
  .train-grid input {
    margin: 0;
  }
  :is(.fs-mode, .filter-row) {
    display: flex;
    gap: 1em;
    align-items: center;
    justify-content: space-between;
  }
  /* label left + vertically centered, radio options stacked on the right (mirrors the
  label-left/inputs-right layout of the require/exclude rows above) */
  .train-grid .fs-mode {
    grid-column: 1 / -1;
    margin-top: 1pt;
  }
  .fs-mode-options {
    display: flex;
    flex-direction: column;
    gap: 2pt;
  }
  .fs-mode-options label {
    display: flex;
    gap: 3pt;
    align-items: center;
  }
  .fs-mode-options label {
    justify-content: end;
  }
  button.clear-filters {
    display: flex;
    gap: 3pt;
    align-items: center;
    background: none;
    padding: 0;
    color: var(--link-color);
  }
  button.preset {
    padding: 1pt 6pt;
    text-align: left;
  }
  .filter-row button.preset {
    flex: 1;
  }
  button.delete-preset {
    background: none;
    padding: 0 2pt;
    opacity: 0.7;
  }
  .dropdown form {
    display: flex;
    gap: 4pt;
    margin-top: 4pt;
  }
  .dropdown form input {
    width: 13em;
    font-size: inherit;
    background: var(--btn-bg);
  }
</style>
