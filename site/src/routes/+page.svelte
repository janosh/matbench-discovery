<script lang="ts">
  import { page } from '$app/state'
  import DiscoverySetToggle from '$lib/DiscoverySetToggle.svelte'
  import MetricsTable from '$lib/table/MetricsTable.svelte'
  import { DISCOVERY_SETS } from '$lib/types'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import GitHubActivityScatter from '$lib/plot/GitHubActivityScatter.svelte'
  import RadarChart from '$lib/plot/RadarChart.svelte'
  import {
    ALL_METRICS,
    DIATOMICS_METRICS,
    DISCOVERY_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    HYPERPARAMS,
    MD_METRICS,
    METADATA_COLS,
    PHONON_METRICS,
  } from '$lib/labels'
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined-scores.svelte'
  import { is_finite_num, metric_value } from '$lib/metrics'
  import { make_table_filters, ACTIVE_MODELS } from '$lib/models.svelte'
  import {
    apply_weights_param,
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    weights_to_param,
  } from '$lib/url-state.svelte'
  import { valid_query_param } from 'svelte-widgets/url-params'
  import type { DiscoverySet, Label, ModelData, SortDir } from '$lib/types'
  import { ButtonGroup, Icon } from 'svelte-widgets'
  import { RSS } from 'svelte-widgets/icons'
  import { onMount } from 'svelte'
  import { slide } from 'svelte/transition'
  import { tooltip } from 'svelte-widgets/attachments'
  import type { Snapshot } from './$types'
  import github_activity_data from './models/mlip-github-activity.json'

  // landing hid TPR; keep its Recall replacement supplementary too
  const supplementary_hidden = new Set(
    [DISCOVERY_METRICS.Recall, DISCOVERY_METRICS.RMSE].map((metric) => metric.key),
  )
  const metadata_keys = new Set([
    ...Object.values(METADATA_COLS).map((col) => col.key),
    HYPERPARAMS.model_params.key,
    HYPERPARAMS.graph_construction_radius.key,
  ])
  const col_presets = {
    Discovery: Object.values(DISCOVERY_METRICS),
    'Geo Opt': Object.values(GEO_OPT_SYMMETRY_METRICS),
    Phonons: Object.values(PHONON_METRICS),
    MD: Object.values(MD_METRICS),
    Diatomics: Object.values(DIATOMICS_METRICS),
  }
  type ColPreset = keyof typeof col_presets
  const default_col_preset: ColPreset = `Discovery`
  const col_preset_names = Object.keys(col_presets) as ColPreset[]
  const preset_expansions: Partial<Record<ColPreset, string>> = {
    'Geo Opt': `Geometry Optimization`,
    MD: `Molecular Dynamics`,
  }
  const preset_primary_metrics: Record<ColPreset, Label> = {
    Discovery: DISCOVERY_METRICS.F1,
    'Geo Opt': ALL_METRICS.RMSD,
    Phonons: PHONON_METRICS.κ_SRME,
    MD: MD_METRICS.md_combined_score,
    Diatomics: DIATOMICS_METRICS.diatomics_combined_score,
  }
  // CPS and one headline metric per task stay visible across presets.
  const headline_metric_keys = new Set(
    [ALL_METRICS.CPS, ...Object.values(preset_primary_metrics)].map(
      (metric) => metric.key,
    ),
  )
  const preset_default_sorts = Object.fromEntries(
    Object.entries(preset_primary_metrics).map(([preset, primary_metric]) => {
      const sort_metric = preset === `Discovery` ? ALL_METRICS.CPS : primary_metric
      return [
        preset,
        {
          column: sort_metric.key,
          dir: sort_metric.better === `lower` ? `asc` : `desc`,
        },
      ]
    }),
  ) as Record<ColPreset, { column: string; dir: SortDir }>
  const filters = make_table_filters()
  const col_preset_options = col_preset_names.map((name) => ({
    value: name,
    label: name,
    tooltip: `Focus the table on ${preset_expansions[name] ?? name} metrics`,
  }))

  let col_preset = $state<ColPreset>(default_col_preset)
  let preset_metric_keys = $derived(
    new Set([...headline_metric_keys, ...col_presets[col_preset].map((col) => col.key)]),
  )
  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
  let sort = $state({ ...preset_default_sorts[default_col_preset] })
  let auto_sort_enabled = $state(true)
  let custom_col_config = $state(false)
  let previous_col_preset: ColPreset = default_col_preset
  const sortable_header_selector = `thead th[role="button"]`
  const column_toggle_input_selector = `.column-menu input[type="checkbox"]`
  const reset_columns_selector = `button[aria-label="Reset all columns to defaults"]`

  $effect(() => {
    if (col_preset === previous_col_preset) return
    previous_col_preset = col_preset
    custom_col_config = false
    if (!auto_sort_enabled) return
    sort = { ...preset_default_sorts[col_preset] }
  })

  function handle_table_event(event: Event) {
    if (event instanceof KeyboardEvent && ![`Enter`, ` `].includes(event.key)) return
    const target = event.target
    if (!(target instanceof Element)) return
    if (target.closest(sortable_header_selector)) auto_sort_enabled = false
    if (target.closest(reset_columns_selector)) custom_col_config = false
    if (target.matches(column_toggle_input_selector)) custom_col_config = true
  }

  const valid_sets = new Set(DISCOVERY_SETS)
  onMount(() => {
    const params = page.url.searchParams
    const next_preset =
      col_preset_names.find((preset) => preset === params.get(`preset`)) ??
      default_col_preset
    const next_sort = sort_from_query(params, preset_default_sorts[next_preset])
    const default_sort = preset_default_sorts[next_preset]
    auto_sort_enabled =
      next_sort.column === default_sort.column && next_sort.dir === default_sort.dir
    sort = next_sort

    discovery_set = valid_query_param(params, `set`, `unique_prototypes`, valid_sets)
    filters.read(params)
    col_preset = next_preset
    previous_col_preset = next_preset
  })

  // Sync table state back to URL query params after the initial URL read (table state
  // is read once in onMount above; weights re-read on every navigation so same-route
  // navs to a weights-less `/` reset them like the MD page does).
  bind_url_params(
    (params) => {
      apply_weights_param(params.get(`weights`), CPS_CONFIG, DEFAULT_CPS_CONFIG)
    },
    () => [
      // omit `preset` for the default and when the user customized
      // columns (a preset no longer describes the visible column set)
      [`preset`, custom_col_config ? default_col_preset : col_preset, default_col_preset],
      [`set`, discovery_set, `unique_prototypes`],
      ...sort_url_entries(sort, preset_default_sorts[col_preset]),
      ...filters.url_entries,
      // custom CPS weights (F1,κ_SRME,RMSD); omitted at defaults
      [`weights`, weights_to_param(CPS_CONFIG, DEFAULT_CPS_CONFIG)],
    ],
  )

  // Each task view includes only models with its headline metric.
  let has_preset_data = $derived((model: ModelData) =>
    is_finite_num(metric_value(model, preset_primary_metrics[col_preset], discovery_set)),
  )
  let in_cohort = $derived(
    (model: ModelData) => has_preset_data(model) && filters.matches(model),
  )

  export const snapshot: Snapshot = {
    capture: () => ({
      discovery_set,
      col_preset,
      custom_col_config,
      sort,
      auto_sort_enabled,
      filters: filters.as_preset,
      show_heatmap: filters.show_heatmap,
    }),
    // Snapshots outlive deploys, so a stored discovery set, column preset or dataset
    // key may since have been renamed or removed. Restoring one would filter every
    // model away or leave a toggle with nothing selected, and the user has no way to
    // see why. So each value is only restored if it still names something real,
    // otherwise the freshly-mounted default stands as if no snapshot existed.
    restore: (values) => {
      custom_col_config = values.custom_col_config ?? custom_col_config
      auto_sort_enabled = values.auto_sort_enabled ?? auto_sort_enabled
      sort = values.sort ?? sort
      if (valid_sets.has(values.discovery_set)) discovery_set = values.discovery_set
      col_preset =
        col_preset_names.find((preset) => preset === values.col_preset) ?? col_preset
      previous_col_preset = col_preset
      // apply() drops unknown dataset keys, targets and openness values itself
      if (values.filters) filters.apply(values.filters)
      filters.show_heatmap = values.show_heatmap ?? filters.show_heatmap
    },
  }
</script>

<!-- MatterViz portals column toggle inputs to document.body, outside the table section. -->
<svelte:document onclickcapture={handle_table_event} />

<h1 id="matbench-discovery">
  <img src="/favicon.svg" alt="Matbench Discovery Logo" width="60px" />
  Matbench Discovery
</h1>

<p class="intro">
  Compare machine-learning models across <a href="/tasks">five materials-science tasks</a
  >.
</p>

<figure id="metrics-table">
  <div class="toggle-row">
    <span>Column presets:</span>
    <ButtonGroup
      bind:selected={col_preset}
      label="Column presets"
      options={col_preset_options}
      tooltip_options={{ placement: `top` }}
    />
  </div>
  <!-- the test-set selector only affects discovery metrics, so only show it in the
  Discovery preset where those columns are visible -->
  {#if col_preset === `Discovery`}
    <div class="toggle-row" in:slide={{ duration: 250 }}>
      <span>Discovery test set:</span>
      <DiscoverySetToggle bind:selected={discovery_set} />
    </div>
  {/if}

  <!-- surface the MD beta warning right at the table when MD columns are shown -->
  {#if col_preset === `MD`}
    <p class="task-note">
      <strong>MD is in beta.</strong> These metrics are preliminary and may change.
      <a href="/tasks/md">About this task →</a>
    </p>
  {/if}

  <section
    class="full-bleed"
    onchangecapture={handle_table_event}
    onkeydowncapture={handle_table_event}
  >
    <MetricsTable
      col_filter={(col) =>
        metadata_keys.has(col.key)
          ? col.visible !== false
          : preset_metric_keys.has(col.key) && !supplementary_hidden.has(col.key)}
      {discovery_set}
      model_filter={has_preset_data}
      bind:sort
      {filters}
    />
  </section>

  <figcaption>
    <div class="table-footer">
      <p>{ACTIVE_MODELS.filter(in_cohort).length} models</p>
      <div style="display: flex; align-items: center; gap: 1em">
        <a href="/contribute">Submit a model</a>
        <a
          href="/rss.xml"
          title="Follow new model submissions in your RSS reader"
          {@attach tooltip()}
        >
          <Icon icon={RSS} /> RSS
        </a>
      </div>
    </div>
    <section id="score-weights" aria-labelledby="score-weights-heading">
      <h3 id="score-weights-heading">Adjust score weights</h3>
      <div class="score-guide">
        <p>
          CPS combines <a href="/tasks/discovery">discovery (F1)</a>,
          <a href="/tasks/geo-opt">geometry optimization (RMSD)</a>, and
          <a href="/tasks/phonons">thermal conductivity (κ<sub>SRME</sub>)</a>. Drag the
          dot to change their importance; scores and rankings update immediately. Custom
          weights are included in the page URL so you can share your view.
        </p>
        <RadarChart size={260} />
      </div>
    </section>
    <details class="page-details" id="table-guide">
      <summary>How to read the table</summary>
      <p>
        Select a column heading to sort. Hover labels for definitions and n/a cells for
        missing-result explanations. Use Compare or double-click model rows to compare
        models side by side.
      </p>
      <p>
        <a href="/data/sets">Training Set</a> counts distinct materials, with relaxation
        frames in parentheses. When only frame counts are available, those are shown
        instead.
        <code>(N=x)</code> beside Params gives the number of estimators in an ensemble.
      </p>
      <p><a href="/tasks">Task definitions and methodology →</a></p>
    </details>
  </figcaption>
</figure>

<section class="plot-section" aria-labelledby="cps-progress-over-time">
  <h2 id="cps-progress-over-time">CPS Progress Over Time</h2>
  <p>Each dot is a model; the dashed line tracks the best score so far.</p>
  <DynamicScatter
    models={ACTIVE_MODELS}
    model_filter={in_cohort}
    {discovery_set}
    x_key={METADATA_COLS.benchmark_added.key}
    y_key={ALL_METRICS.CPS.key}
    show_pareto_frontier
  />
</section>

<section class="plot-section" aria-labelledby="github-activity">
  <h2 id="github-activity">GitHub Activity</h2>
  <p>
    Explore <a href="/models">model communities</a>. Larger dots mean more contributors;
    color shows recent commits.
  </p>
  <GitHubActivityScatter github_data={github_activity_data} />
</section>

<section id="about-benchmark" aria-labelledby="about-matbench-discovery">
  <h2 id="about-matbench-discovery">About Matbench Discovery</h2>
  <p>
    This benchmark compares accuracy, robustness, and computational cost across
    <a href="/tasks/discovery">crystal discovery</a>,
    <a href="/tasks/geo-opt">geometry optimization</a>,
    <a href="/tasks/phonons">phonons</a>,
    <a href="/tasks/md">molecular dynamics</a>, and
    <a href="/tasks/diatomics">diatomics</a>. Rankings help you explore trade-offs; they
    are not a complete assessment or an endorsement of a model.
  </p>
  <p>
    Crystal stability is evaluated against a
    <a href="/tasks/discovery#convex-hull-construction-in-matbench-discovery"
      >convex hull</a
    >
    built from
    <a
      href="https://docs.materialsproject.org/methodology/materials-methodology/calculation-details"
      >DFT</a
    >
    reference energies.
  </p>
  <p>
    Riebesell, J., Goodall, R.E.A., Benner, P. et al.
    <a href="https://doi.org/10.1038/s42256-025-01055-1">
      A framework to evaluate machine learning crystal stability predictions.</a
    >
    <i>Nature Machine Intelligence</i> 7, 836–847 (2025).
  </p>
  <p>
    <a href="/tasks">Explore the tasks</a> ·
    <a href="/contribute">Contribute predictions</a>
  </p>
</section>

<style>
  h1 {
    margin-block: -1.2em 0.35em;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 7pt;
  }
  h1 img {
    filter: brightness(0.8);
  }
  :root[data-theme='light'] h1 img {
    filter: brightness(0.2);
  }
  figure {
    margin: 0;
    display: grid;
    gap: 1ex;
  }
  :is(.toggle-row, .table-footer) {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
  }
  .toggle-row {
    gap: 8pt;
    font-size: smaller;
  }
  .intro {
    text-align: center;
    color: var(--text-muted);
    margin: 0 0 1.8em;
  }
  figcaption {
    font-size: 0.9em;
  }
  .table-footer {
    justify-content: space-between;
    gap: 0 1em;
    > p {
      flex: 1 1 25em;
      color: var(--text-muted);
      margin-block: 0.75em;
    }
  }
  #score-weights,
  .page-details {
    border-top: 1px solid var(--border);
    padding-block: 0.65em;
    > summary {
      color: var(--link-color);
      width: fit-content;
    }
  }
  .score-guide {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
    gap: 1em;
    > p {
      flex: 1 1 18em;
    }
  }
  .plot-section {
    margin-block-start: 2.5em;
    > h2 {
      margin-block-end: 0.4em;
    }
    > p {
      text-align: center;
      color: var(--text-muted);
      margin-block-start: 0;
    }
  }
</style>
