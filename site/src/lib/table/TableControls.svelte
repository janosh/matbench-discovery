<script lang="ts">
  import type { Label } from '$lib'
  import { openness_tooltips } from '$lib/metrics'
  import { make_table_filters } from '$lib/models.svelte'
  import { OPENNESS_OPTIONS, type UrlTableFilters } from '$lib/url-state.svelte'
  import { Icon } from 'matterviz'
  import { type Label as MvLabel, ToggleMenu } from 'matterviz/table'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    show_energy_only = $bindable(false),
    columns = $bindable([]),
    filters = make_table_filters(),
    show_selected_only = $bindable(false),
    selected_count = 0,
    show_energy_only_toggle = false,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    columns?: Label[]
    show_energy_only?: boolean
    filters?: UrlTableFilters
    show_selected_only?: boolean
    selected_count?: number
    show_energy_only_toggle?: boolean
  } = $props()

  const close_on_outside_click = click_outside({
    callback: (node) => ((node as HTMLDetailsElement).open = false),
  })

  const train_modes = [
    { mode: `require`, label: `only`, title: `Only show models trained on` },
    { mode: `exclude`, label: `not`, title: `Only show models NOT trained on` },
  ] as const
  const n_train = $derived(Object.keys(filters.training).length)
  const n_openness = $derived(filters.openness.length)
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
      {@attach tooltip()}
    >
      Training data{n_train ? ` (${n_train})` : ``}
    </summary>
    <div class="dropdown">
      <span class="hint">
        <em>only</em> = only models trained on this dataset (combined constraints must all
        hold), <em>not</em> = only models <em>not</em> trained on it
      </span>
      {#each filters.training_sets as dataset_key (dataset_key)}
        <div class="filter-row">
          <span>{dataset_key}</span>
          {#each train_modes as { mode, label, title } (mode)}
            <label title="{title} {dataset_key}" {@attach tooltip()}>
              <input
                type="checkbox"
                checked={filters.training[dataset_key] === mode}
                onchange={() => filters.set_training(dataset_key, mode)}
              />
              {label}
            </label>
          {/each}
        </div>
      {/each}
    </div>
  </details>

  <details class="filter-menu" {@attach close_on_outside_click}>
    <summary
      title="Filter models by whether their source code and training data are open"
      {@attach tooltip()}
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
          {openness}
        </label>
      {/each}
    </div>
  </details>

  {#if filters.n_active > 0}
    <button
      class="clear-filters"
      onclick={() => filters.clear()}
      title="Reset training-data and openness filters"
      {@attach tooltip()}
    >
      <Icon icon="Cross" /> clear filters
    </button>
  {/if}

  {#if show_energy_only_toggle}
    <label>
      <input type="checkbox" bind:checked={show_energy_only} />
      Energy-only models
      <span
        title="Include models that only predict energy (no forces or stress)"
        {@attach tooltip()}
      >
        <Icon icon="Info" />
      </span>
    </label>
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
  div.table-controls {
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
    gap: 3pt;
    min-width: max-content;
    margin-top: 4px;
    padding: 6pt 8pt;
    background: var(--page-bg);
    border: 1px solid var(--border);
    border-radius: 5px;
    box-shadow: 0 0 10px var(--shadow);
  }
  details.filter-menu .hint {
    max-width: 28em;
    font-size: 0.85em;
    opacity: 0.75;
    text-wrap: balance;
    margin-bottom: 3pt;
  }
  .filter-row {
    display: flex;
    gap: 1em;
    align-items: center;
    justify-content: space-between;
  }
  .filter-row > span:first-child {
    margin-right: auto;
  }
  .filter-row label {
    display: flex;
    gap: 3pt;
    align-items: center;
  }
  button.clear-filters {
    display: flex;
    gap: 3pt;
    align-items: center;
    background: none;
    padding: 0;
    color: var(--link-color);
  }
</style>
