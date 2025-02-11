<script lang="ts">
  import {
    MetricScatter,
    MetricsTable,
    MODEL_METADATA,
    TableColumnToggleMenu,
  } from '$lib'
  import { METADATA_COLS, METRICS_COLS } from '$lib/metrics'
  import type { DiscoverySet } from '$lib/types'
  import Icon from '@iconify/svelte'
  import { Tooltip } from 'svelte-zoo'
  import { click_outside } from 'svelte-zoo/actions'

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
    (md) => md.metrics?.discovery?.[discovery_set]?.F1 != null,
  )

  // Get array of hidden columns
  $: hide_cols = Object.entries(visible_cols)
    .filter(([_, visible]) => !visible)
    .map(([col]) => col)
</script>

<h1>Crystal Stability Prediction Metrics</h1>

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

  <MetricsTable {hide_cols} {discovery_set} style="width: 100%;" />

  <div class="table-controls">
    <TableColumnToggleMenu bind:visible_cols bind:column_panel_open />
  </div>

  <h3>F1 classification score of models over time</h3>
  <p>
    The F1 score is the harmonic mean of precision and recall. It is a measure of the
    model's ability to correctly identify hypothetical crystals in the WBM test set as
    lying on or below the Materials Project convex hull.
  </p>

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
  h3 {
    text-align: center;
  }
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
