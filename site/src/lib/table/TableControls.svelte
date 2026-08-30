<script lang="ts">
  import type { TableLabel } from '$lib'
  import {
    BUILTIN_PRESETS,
    delete_user_preset,
    save_user_preset,
    user_presets,
  } from '$lib/filter-presets.svelte'
  import { openness_tooltips } from '$lib/metrics'
  import { comparison } from '$lib/model-comparison.svelte'
  import { make_table_filters, ACTIVE_MODELS } from '$lib/models.svelte'
  import {
    DEFAULT_TARGETS_PARAM,
    FS_MODES,
    OPENNESS_OPTIONS,
    parse_targets,
    TARGET_OUTPUTS,
    TRAIN_FILTER_MODES,
  } from '$lib/url-state.svelte'
  import type { TargetOutput, UrlTableFilters } from '$lib/url-state.svelte'
  import { Icon, Sheet } from 'svelte-widgets'
  import { Cross, Filter, Scale } from 'svelte-widgets/icons'
  import { ToggleMenu } from 'matterviz/table'
  import type { Snippet } from 'svelte'
  import { click_outside, tooltip } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    columns = $bindable([]),
    filters = make_table_filters(),
    show_selected_only = $bindable(false),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    columns?: TableLabel[]
    filters?: UrlTableFilters
    show_selected_only?: boolean
  } = $props()
  let selected_count = $derived(comparison.keys.size)

  const close_on_outside_click = click_outside({
    callback: (node) => node.removeAttribute(`open`),
  })

  // summaries open their pane directly below, so tooltips go above the button
  const summary_tooltip = tooltip({ placement: `top` })
  // Sheet open state is CSS-driven for layout; kept in sync so a resize from mobile
  // to desktop can dismiss an open dialog (display:none alone can leave it top-layered)
  let filter_sheet_open = $state(false)

  const target_outputs = Object.entries(TARGET_OUTPUTS) as [TargetOutput, string][]
  // static per-category model tallies shown in the filter panels. Semantics mirror
  // UrlTableFilters.matches: targets use the same parser, and fs_mode `any` counts
  // every model.
  const openness_counts: Record<string, number> = {}
  const training_counts: Record<string, number> = {}
  const target_counts: Record<string, number> = {}
  const fs_mode_counts: Record<string, number> = { any: ACTIVE_MODELS.length }
  for (const model of ACTIVE_MODELS) {
    const openness = model.openness
    openness_counts[openness] = (openness_counts[openness] ?? 0) + 1
    for (const dataset of model.training_sets) {
      training_counts[dataset] = (training_counts[dataset] ?? 0) + 1
    }
    const { outputs, fs_mode } = parse_targets(model.targets)
    for (const output of outputs) {
      target_counts[output] = (target_counts[output] ?? 0) + 1
    }
    if (fs_mode) fs_mode_counts[fs_mode] = (fs_mode_counts[fs_mode] ?? 0) + 1
  }
  const n_train = $derived(Object.keys(filters.training).length)
  const training_sets_by_model_count = $derived(
    filters.training_sets.toSorted(
      (dataset_left, dataset_right) =>
        (training_counts[dataset_right] ?? 0) - (training_counts[dataset_left] ?? 0),
    ),
  )
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
  const close_sheet_on_desktop = (): void => {
    if (globalThis.innerWidth > 600) filter_sheet_open = false
  }
</script>

<svelte:window onresize={close_sheet_on_desktop} />

{#snippet filter_mode_headers()}
  {#each TRAIN_FILTER_MODES as mode, mode_idx (mode)}
    <span class="col-head" style:grid-column={mode_idx + 2}>{mode}</span>
  {/each}
{/snippet}

{#snippet training_filters()}
  <div class="filter-content train-grid">
    <span class="hint">
      <em>require</em> = model's training set must include this dataset,
      <em>exclude</em> = hide models trained on it
    </span>
    {@render filter_mode_headers()}
    {#each training_sets_by_model_count as dataset_key (dataset_key)}
      <span>{dataset_key} ({training_counts[dataset_key] ?? 0})</span>
      {#each TRAIN_FILTER_MODES as mode (mode)}
        <input
          type="checkbox"
          aria-label="{mode} {dataset_key}"
          checked={filters.training[dataset_key] === mode}
          onchange={() => filters.set_training(dataset_key, mode)}
        />
      {/each}
    {/each}
  </div>
{/snippet}

{#snippet openness_filters()}
  <div class="filter-content">
    {#each OPENNESS_OPTIONS as openness (openness)}
      <label class="filter-row" title={openness_tooltips[openness]} {@attach tooltip()}>
        <input
          type="checkbox"
          checked={filters.openness.includes(openness)}
          onchange={(event) => {
            filters.toggle_openness(openness)
            // toggle can refuse to clear the last option; restore native checked
            event.currentTarget.checked = filters.openness.includes(openness)
          }}
        />
        {openness} ({openness_counts[openness] ?? 0})
      </label>
    {/each}
  </div>
{/snippet}

{#snippet target_filters()}
  <div class="filter-content train-grid">
    <span class="hint">
      Every model predicts energy (E). <em>require</em>/<em>exclude</em> filter by the other
      predicted outputs; forces are required by default (hides energy-only models)
    </span>
    {@render filter_mode_headers()}
    {#each target_outputs as [key, label] (key)}
      <span>{label} ({key}) ({target_counts[key] ?? 0})</span>
      {#each TRAIN_FILTER_MODES as mode (mode)}
        <input
          type="checkbox"
          aria-label="{mode} {label}"
          checked={filters.targets[key] === mode}
          onchange={() => filters.set_target(key, mode)}
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
{/snippet}

{#snippet preset_filters()}
  <div class="filter-content">
    {#each Object.entries( { ...BUILTIN_PRESETS, ...user_presets } ) as [name, preset] (name)}
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
            <Icon icon={Cross} />
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
{/snippet}

{#snippet filter_section(mobile: boolean, label: string, title: string, content: Snippet)}
  {#if mobile}
    <section class="sheet-section">
      <h3>{label}</h3>
      {@render content()}
    </section>
  {:else}
    <details class="filter-menu" {@attach close_on_outside_click}>
      <summary {title} {@attach summary_tooltip}>{label}</summary>
      {@render content()}
    </details>
  {/if}
{/snippet}

{#snippet filter_sections(mobile: boolean)}
  {@render filter_section(
    mobile,
    `Training data${n_train ? ` (${n_train})` : ``}`,
    `Filter models by the datasets they were trained on`,
    training_filters,
  )}
  {@render filter_section(
    mobile,
    `Openness${filters.openness.length < OPENNESS_OPTIONS.length ? ` (${filters.openness.length}/${OPENNESS_OPTIONS.length})` : ``}`,
    `Filter models by whether their source code and training data are open`,
    openness_filters,
  )}
  {@render filter_section(
    mobile,
    `Targets${targets_badge}`,
    `Filter models by which quantities they predict and how forces/stress are computed`,
    target_filters,
  )}
  {@render filter_section(
    mobile,
    `Presets`,
    `Apply a saved filter combination or save the current one`,
    preset_filters,
  )}
{/snippet}

<div class="table-controls" {...rest}>
  <button
    class="compare"
    onclick={() => (comparison.open = true)}
    title={selected_count
      ? `Open the side-by-side comparison of the ${selected_count} selected model${selected_count === 1 ? `` : `s`}`
      : `Compare models side by side. Select them by hovering a row and clicking its ⚖ button, right-clicking a row, or double-clicking it`}
    {@attach tooltip()}
  >
    <Icon icon={Scale} /> Compare{selected_count ? ` (${selected_count})` : ``}
  </button>
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

  <!-- Both trees always mount; CSS media queries toggle visibility. A JS MediaQuery
  branch (SSR fallback false → desktop <details>, mobile client → <Sheet>) hydrated
  mismatched markup. -->
  <div class="mobile-filters">
    <Sheet
      bind:open={filter_sheet_open}
      aria-label="Model filters"
      style="--sheet-size: min(26rem, 100vw); --sheet-content-padding: 0.75rem;"
    >
      {#snippet trigger(trigger_props)}
        <button class="filter-sheet-trigger" type="button" {...trigger_props}>
          <Icon icon={Filter} /> Filters{filters.n_active ? ` (${filters.n_active})` : ``}
        </button>
      {/snippet}
      {#snippet header({ close })}
        <div class="filter-sheet-header">
          <strong>Model filters</strong>
          <button aria-label="Close model filters" onclick={close} type="button">
            <Icon icon={Cross} />
          </button>
        </div>
      {/snippet}
      {@render filter_sections(true)}
    </Sheet>
  </div>
  <div class="desktop-filters">
    {@render filter_sections(false)}
  </div>

  {#if filters.n_active > 0}
    <button
      class="clear-filters"
      onclick={() => filters.clear()}
      title="Reset training-data, openness and target filters"
      {@attach tooltip()}
    >
      <Icon icon={Cross} /> clear filters
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

  <ToggleMenu bind:columns />
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
  .mobile-filters {
    display: none;
  }
  .desktop-filters {
    display: contents;
  }
  details.filter-menu {
    position: relative;
    summary {
      list-style: none;
      padding: 1pt 6pt;
      border-radius: 4px;
      background: var(--btn-bg);
    }
    &[open] summary {
      background: color-mix(in srgb, var(--link-color) 25%, transparent);
    }
    .filter-content {
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
  }
  .hint {
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
    .hint {
      grid-column: 1 / -1;
      /* contribute zero to grid column sizing (the dataset rows set the panel width),
      then stretch to whatever width they produce and wrap */
      width: 0;
      min-width: 100%;
      text-wrap: wrap;
    }
    :is(.col-head, input) {
      justify-self: center;
    }
    .col-head {
      font-style: italic;
      opacity: 0.75;
    }
    input {
      margin: 0;
    }
    .fs-mode {
      grid-column: 1 / -1;
      margin-top: 1pt;
    }
  }
  :is(.fs-mode, .filter-row) {
    display: flex;
    gap: 1em;
    align-items: center;
    justify-content: space-between;
  }
  /* label left + vertically centered, radio options stacked on the right (mirrors the
  label-left/inputs-right layout of the require/exclude rows above) */
  .fs-mode-options {
    display: flex;
    flex-direction: column;
    gap: 2pt;
    label {
      display: flex;
      gap: 3pt;
      align-items: center;
      justify-content: end;
    }
  }
  :is(button.clear-filters, button.filter-sheet-trigger, button.compare) {
    display: inline-flex;
    gap: 3pt;
    align-items: center;
    color: var(--link-color);
  }
  button.clear-filters {
    background: none;
    padding: 0;
  }
  /* the one action here that isn't a filter: filled so it reads as primary, and set off
     from the filter cluster by a divider */
  button.compare {
    position: relative;
    padding: 1pt 8pt;
    border-radius: 1em;
    background: var(--link-color);
    color: var(--page-bg);
    font-weight: 600;
    &::after {
      content: '';
      position: absolute;
      right: -6pt; /* centered in the 12pt column gap */
      height: 1.2em;
      border-left: 1px solid var(--border);
    }
  }
  .filter-row button.preset {
    flex: 1;
    padding: 1pt 6pt;
    text-align: left;
  }
  button.delete-preset {
    background: none;
    padding: 0 2pt;
    opacity: 0.7;
  }
  .filter-content form {
    display: flex;
    gap: 4pt;
    margin-top: 4pt;
    input {
      width: 13em;
      font-size: inherit;
      background: var(--btn-bg);
    }
  }
  .filter-sheet-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    button {
      display: grid;
      place-items: center;
      padding: 3pt;
      background: transparent;
    }
  }
  .sheet-section {
    padding-block: 0 0.9rem;
    border-bottom: 1px solid var(--border);
    &:last-child {
      border-bottom: 0;
    }
    h3 {
      margin-block: 0 0.45rem;
      font-size: 1rem;
    }
    .filter-content {
      display: grid;
      &:not(.train-grid) {
        gap: 4pt;
      }
    }
  }
  @media (max-width: 600px) {
    .mobile-filters {
      display: contents;
    }
    .desktop-filters {
      display: none;
    }
  }
</style>
