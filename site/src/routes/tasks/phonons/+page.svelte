<script lang="ts">
  import { MetricsTable, model_is_compliant, MODEL_METADATA } from '$lib'
  import { METADATA_COLS, METRICS_COLS } from '$lib/metrics'
  import MetricScatter from '$lib/MetricScatter.svelte'
  import ColumnToggle from '$site/src/lib/TableColumnToggleMenu.svelte'
  import Icon from '@iconify/svelte'
  import { Toggle, Tooltip } from 'svelte-zoo'

  let show_non_compliant: boolean = false
  let show_energy_only: boolean = false

  // Default column visibility
  let visible_cols: Record<string, boolean> = {
    ...Object.fromEntries([...METRICS_COLS].map((col) => [col.label, false])),
    ...Object.fromEntries([...METADATA_COLS].map((col) => [col.label, true])),
    'κ<sub>SRME</sub>': true,
  }

  let kappa_tooltip_point: { x: number; y: number } | null = null
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

<h1>MLFF Phonon Modeling Metrics</h1>

<figure>
  <MetricsTable
    {show_non_compliant}
    {hide_cols}
    {show_energy_only}
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
    metric="phonons.kappa_103.κ_SRME"
    y_label="kappa SRME (lower better)"
    bind:tooltip_point={kappa_tooltip_point}
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
  div.table-controls {
    display: flex;
    flex-wrap: wrap;
    gap: 5pt;
    place-content: center;
    margin: 3pt auto;
  }
</style>
