<script lang="ts">
  import {
    MetricScatter,
    MetricsTable,
    MODEL_METADATA,
    TableColumnToggleMenu,
  } from '$lib'
  import { DISCOVERY_METRICS, DISCOVERY_SET_LABELS, METADATA_COLS } from '$lib/metrics'
  import type { DiscoverySet } from '$lib/types'
  import 'iconify-icon'
  import { Tooltip } from 'svelte-zoo'
  import { click_outside } from 'svelte-zoo/actions'

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    ...Object.fromEntries(
      [...DISCOVERY_METRICS, ...METADATA_COLS].map((col) => [col.label, true]),
    ),
    'κ<sub>SRME</sub>': false,
  })

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
  let f1_tooltip_point: { x: number; y: number } | null = $state(null)
  let hovered = $state(false)
  let column_panel_open: boolean = $state(false)

  let filtered_models = $derived(
    Object.values(MODEL_METADATA).filter(
      (md) => md.metrics?.discovery?.[discovery_set]?.F1 != null,
    ),
  )
</script>

<h1>Crystal Stability Prediction Metrics</h1>

<figure>
  <div
    class="discovery-set-toggle"
    use:click_outside={{ callback: () => (column_panel_open = false) }}
  >
    {#each Object.entries(DISCOVERY_SET_LABELS) as [key, { title, tooltip, link }]}
      <Tooltip text={tooltip} tip_style="z-index: 2; font-size: 0.8em;">
        <button
          class:active={discovery_set === key}
          onclick={() => (discovery_set = key)}
        >
          {title}
          {#if link}
            <a
              href={link}
              aria-label="Open in new tab"
              target="_blank"
              rel="noopener noreferrer"
            >
              <iconify-icon icon="octicon:info" inline></iconify-icon>
            </a>
          {/if}
        </button>
      </Tooltip>
    {/each}
  </div>

  <MetricsTable
    col_filter={(col) => visible_cols[col.label] ?? true}
    {discovery_set}
    style="width: 100%;"
  />

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
