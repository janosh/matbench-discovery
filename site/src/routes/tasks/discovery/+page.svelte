<script lang="ts">
  import { MetricsTable, model_is_compliant, MODEL_METADATA } from '$lib'
  import { METADATA_COLS, METRICS_COLS } from '$lib/metrics'
  import MetricScatter from '$lib/MetricScatter.svelte'
  import type { DiscoverySet } from '$lib/types'
  import ColumnToggle from '$site/src/lib/TableColumnToggleMenu.svelte'
  import Icon from '@iconify/svelte'
  import { Toggle, Tooltip } from 'svelte-zoo'
  import { click_outside } from 'svelte-zoo/actions'

  let show_non_compliant: boolean = false
  let show_energy_only: boolean = false

  // Default column visibility
  let visible_cols: Record<string, boolean> = {
    ...Object.fromEntries(
      [...METRICS_COLS, ...METADATA_COLS].map((col) => [col.label, true]),
    ),
    'κ<sub>SRME</sub>': false,
  }

  const discovery_set_labels: Record<
    DiscoverySet,
    { title: string; tooltip: string; link?: string }
  > = {
    full_test_set: {
      title: `Full Test Set`,
      tooltip: `Metrics computed on the full test set including duplicate structure prototypes`,
    },
    unique_prototypes: {
      title: `Unique Prototypes`,
      tooltip: `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
      link: `https://github.com/janosh/matbench-discovery/blob/fd1dda6c/data/wbm/compile_wbm_test_set.py#L632-L705`,
    },
    most_stable_10k: {
      title: `10k Most Stable`,
      tooltip: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
    },
  }
  let discovery_set: DiscoverySet = `unique_prototypes`

  let f1_tooltip_point: { x: number; y: number } | null = null
  let hovered = false
  let column_panel_open: boolean = false

  $: filtered_models = Object.values(MODEL_METADATA).filter(
    (md) => show_non_compliant || model_is_compliant(md),
  )

  // Get array of hidden columns
  $: hide_cols = Object.entries(visible_cols)
    .filter(([_, visible]) => !visible)
    .map(([col]) => col)
</script>

<h1>MLFF Materials Discovery Metrics</h1>

<figure>
  <div
    class="discovery-set-toggle"
    use:click_outside={{ callback: () => (column_panel_open = false) }}
  >
    {#each Object.entries(discovery_set_labels) as [key, { title, tooltip, link }]}
      <Tooltip text={tooltip} tip_style="z-index: 2; font-size: 0.8em;">
        <button
          class:active={discovery_set === key}
          on:click={() => (discovery_set = key)}
        >
          {title}
          {#if link}
            <a href={link} target="_blank">
              <Icon icon="octicon:info" inline />
            </a>
          {/if}
        </button>
      </Tooltip>
    {/each}
  </div>

  <MetricsTable
    {show_non_compliant}
    {hide_cols}
    {show_energy_only}
    {discovery_set}
    style="width: 100%;"
  />

  <div class="table-controls">
    <Toggle bind:checked={show_non_compliant} style="gap: 3pt;">
      Show non-compliant models <Tooltip max_width="20em">
        <span slot="tip">
          Models can be non-compliant for multiple reasons<br />
          - closed source (model implementation and/or train/test code)<br />
          - closed weights<br />
          - trained on more than the permissible training set (<a
            href="https://docs.materialsproject.org/changes/database-versions#v2022.10.28"
            >MP v2022.10.28 release</a
          >)<br />
          We still show these models behind a toggle as we expect them<br /> to nonetheless
          provide helpful signals for developing future models.
        </span>
        <Icon icon="octicon:info-16" inline style="padding: 0 3pt;" />
      </Tooltip>&ensp;</Toggle
    >
    <Toggle bind:checked={show_energy_only} style="gap: 3pt;">
      Show energy-only models <Tooltip max_width="12em">
        <span slot="tip">
          Models that only predict energy (E) perform worse<br /> and can't be evaluated
          on force-modeling tasks such as κ<sub>SRME</sub>
        </span>
        <Icon icon="octicon:info-16" inline style="padding: 0 3pt;" />
      </Tooltip>&ensp;</Toggle
    >

    <ColumnToggle bind:visible_cols bind:column_panel_open />
  </div>

  <MetricScatter
    models={filtered_models}
    metric={`discovery.${discovery_set}.F1`}
    y_label="F1 Score (higher better)"
    bind:tooltip_point={f1_tooltip_point}
    bind:hovered
    style="margin: 2em 0;"
  />
</figure>

<style>
  figure {
    margin: 0;
    display: grid;
    gap: 1ex;
  }
  .discovery-set-toggle {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 5pt;
    margin-bottom: 5pt;
  }
  .discovery-set-toggle button {
    padding: 4px 8px;
    border: 1px solid rgba(255, 255, 255, 0.05);
    background: transparent;
  }
  .discovery-set-toggle button:hover {
    background: rgba(255, 255, 255, 0.05);
  }
  .discovery-set-toggle button.active {
    background: rgba(255, 255, 255, 0.1);
  }
  div.table-controls {
    display: flex;
    flex-wrap: wrap;
    gap: 5pt;
    place-content: center;
    margin: 3pt auto;
  }
</style>
