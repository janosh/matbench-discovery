<script lang="ts">
  import { page } from '$app/state'
  import { DATASETS, DISCOVERY_SETS, MetricsTable, SelectToggle } from '$lib'
  import { DynamicScatter, GitHubActivityScatter, RadarChart } from '$lib/plot'
  import {
    ALL_METRICS,
    DIATOMICS_METRICS,
    DISCOVERY_METRICS,
    discovery_set_toggle_options,
    GEO_OPT_SYMMETRY_METRICS,
    HYPERPARAMS,
    MD_METRICS,
    METADATA_COLS,
  } from '$lib/labels'
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined-scores.svelte'
  import { get_nested_number, is_finite_num, label_data_path } from '$lib/metrics'
  import { make_table_filters, ACTIVE_MODELS } from '$lib/models.svelte'
  import {
    apply_weights_param,
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
    weights_to_param,
  } from '$lib/url-state.svelte'
  import {
    generate_csv,
    generate_excel,
    generate_png,
    generate_svg,
    handle_export,
  } from '$lib/table-export'
  import type { DiscoverySet, Label, ModelData, SortDir } from '$lib/types'
  import Readme from '$root/readme.md'
  import MdNote from '$routes/tasks/md/md-note.md'
  import KappaNote from '$routes/tasks/phonons/kappa-note.md'
  import { format_num, Icon } from 'matterviz'
  import { onMount } from 'svelte'
  import { slide } from 'svelte/transition'
  import { tooltip } from 'svelte-multiselect/attachments'
  import type { Snapshot } from './$types'
  import github_activity_data from './models/mlip-github-activity.json'

  const n_wbm_stable_uniq_protos = 32_942
  const n_wbm_uniq_protos = DATASETS.WBM.n_materials

  const supplementary_hidden = new Set(
    [DISCOVERY_METRICS.TPR, DISCOVERY_METRICS.TNR, DISCOVERY_METRICS.RMSE].map(
      (metric) => metric.key,
    ),
  )
  const metadata_keys = new Set([
    ...Object.values(METADATA_COLS).map((col) => col.key),
    HYPERPARAMS.model_params.key,
    HYPERPARAMS.graph_construction_radius.key,
  ])
  const col_presets = {
    Discovery: Object.values(DISCOVERY_METRICS),
    'Geo Opt': Object.values(GEO_OPT_SYMMETRY_METRICS),
    Phonons: [ALL_METRICS.κ_SRE],
    MD: Object.values(MD_METRICS),
    Diatomics: Object.values(DIATOMICS_METRICS),
  }
  type ColPreset = keyof typeof col_presets
  const default_col_preset: ColPreset = `Discovery`
  const col_preset_names = Object.keys(col_presets) as ColPreset[]
  const preset_primary_metrics: Record<ColPreset, Label> = {
    Discovery: DISCOVERY_METRICS.F1,
    'Geo Opt': ALL_METRICS.RMSD,
    Phonons: ALL_METRICS.κ_SRME,
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
    Object.entries(preset_primary_metrics).map(([preset, metric]) => [
      preset,
      { column: metric.key, dir: metric.better === `lower` ? `asc` : `desc` },
    ]),
  ) as Record<ColPreset, { column: string; dir: SortDir }>
  const filters = make_table_filters()
  const col_preset_options = col_preset_names.map((name) => ({
    value: name,
    label: name,
    tooltip: `Focus the table on ${name} metrics`,
  }))

  let export_error: string | null = $state(null)

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
  const column_toggles_selector = `details.column-toggles`
  const reset_columns_selector = `${column_toggles_selector} button[aria-label="Reset all columns to defaults"]`

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
    if (event.type === `change` && target.matches(`${column_toggles_selector} input`)) {
      custom_col_config = true
    }
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

  let export_state = $derived({ export_error, discovery_set })

  const preset_metric_value = (
    model: ModelData,
    preset: ColPreset,
  ): number | undefined => {
    const discovery = model.metrics?.discovery
    const value =
      preset === `Discovery`
        ? discovery?.[discovery_set]?.F1
        : get_nested_number(model, label_data_path(preset_primary_metrics[preset]))
    return is_finite_num(value) ? value : undefined
  }

  // Each task view includes only models with its headline metric.
  let has_preset_data = $derived(
    (model: ModelData) => preset_metric_value(model, col_preset) !== undefined,
  )
  let in_cohort = $derived(
    (model: ModelData) => has_preset_data(model) && filters.matches(model),
  )
  let primary_metric = $derived(preset_primary_metrics[col_preset])
  let best_entry = $derived.by(() => {
    const entries = ACTIVE_MODELS.filter(in_cohort).flatMap((model) => {
      const value = preset_metric_value(model, col_preset)
      return value === undefined ? [] : [{ model, value }]
    })
    const sort_factor = primary_metric.better === `lower` ? 1 : -1
    return entries.toSorted(
      (entry_1, entry_2) => sort_factor * (entry_1.value - entry_2.value),
    )[0]
  })

  export const snapshot: Snapshot = {
    capture: () => ({
      discovery_set,
      col_preset,
      custom_col_config,
      sort,
      auto_sort_enabled,
      training_filter: { ...filters.training },
      openness: [...filters.openness],
      targets: { ...filters.targets },
      fs_mode: filters.fs_mode,
      show_heatmap: filters.show_heatmap,
    }),
    restore: (values) => {
      custom_col_config = values.custom_col_config ?? custom_col_config
      auto_sort_enabled = values.auto_sort_enabled ?? auto_sort_enabled
      sort = values.sort ?? sort
      discovery_set = values.discovery_set ?? discovery_set
      col_preset = values.col_preset ?? col_preset
      previous_col_preset = col_preset
      filters.training = values.training_filter ?? filters.training
      filters.openness = values.openness ?? filters.openness
      filters.targets = values.targets ?? filters.targets
      filters.fs_mode = values.fs_mode ?? filters.fs_mode
      filters.show_heatmap = values.show_heatmap ?? filters.show_heatmap
    },
  }
</script>

<h1>
  <img src="/favicon.svg" alt="Matbench Discovery Logo" width="60px" />
  Matbench Discovery
</h1>

<figure style="margin-top: 3em" id="metrics-table">
  <div class="toggle-row">
    <span>Column presets:</span>
    <SelectToggle
      bind:selected={col_preset}
      options={col_preset_options}
      tooltip_placement="top"
    />
  </div>
  <!-- the test-set selector only affects discovery metrics, so only show it in the
  Discovery preset where those columns are visible -->
  {#if col_preset === `Discovery`}
    <div class="toggle-row" in:slide={{ duration: 250 }}>
      <span>Discovery test set:</span>
      <SelectToggle
        bind:selected={discovery_set}
        options={discovery_set_toggle_options}
      />
    </div>
  {/if}

  <!-- surface the MD beta warning right at the table when MD columns are shown -->
  {#if col_preset === `MD`}
    <MdNote />
  {/if}

  <section
    class="full-bleed"
    onclickcapture={handle_table_event}
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

  {#if export_error}
    <div class="export-error">
      {export_error}
    </div>
  {/if}

  <div class="downloads">
    Download table as
    {#each [[`SVG`, generate_svg], [`PNG`, generate_png], [`CSV`, generate_csv], [`Excel`, generate_excel]] as const as [label, generate_fn] (label)}
      <button
        class="download-btn"
        onclick={handle_export(generate_fn, label, export_state)}
      >
        {label}
      </button>
    {/each}
    &emsp;Subscribe via
    <a
      href="/rss.xml"
      class="download-btn"
      title="Be notified of new model submissions through an RSS reader"
      {@attach tooltip()}
    >
      <Icon icon="RSS" /> RSS
    </a>
  </div>

  <figcaption class="caption-radar-container">
    <div
      style="background-color: var(--card-bg); padding: 0.2em 0.5em; border-radius: 4px"
    >
      The <strong>CPS</strong> (Combined Performance Score) is a metric that weights
      discovery performance (F1), geometry optimization quality (RMSD), and thermal
      conductivity prediction accuracy (κ<sub>SRME</sub>). Use the radar chart to adjust
      the importance of each component.
      <br /><br />
      The training set column shows the number of materials used to train the model. For models
      trained on DFT relaxations, we show the number of distinct frames in parentheses. In cases
      where only the number of frames is known, we report the number of frames as the training
      set size. <code>(N=x)</code> in the Model Params column shows the number of
      estimators if an ensemble was used.
      {#if col_preset === `Discovery`}
        DAF = Discovery Acceleration Factor measures how many more stable materials a
        model finds compared to random selection from the test set. The unique structure
        prototypes in the WBM test set
        {#if n_wbm_uniq_protos}
          have a
          <code>{format_num(n_wbm_stable_uniq_protos / n_wbm_uniq_protos, `.1%`)}</code>
          rate of stable crystals, meaning the max possible DAF is
          <code>
            ({format_num(n_wbm_stable_uniq_protos)} / {format_num(n_wbm_uniq_protos)})^−1
            ≈ {format_num(n_wbm_uniq_protos / n_wbm_stable_uniq_protos)}
          </code>.
        {:else}
          have an unknown rate of stable crystals (WBM n_materials unavailable).
        {/if}
      {/if}
    </div>
    <!-- CPS weight controls -->
    <RadarChart size={260} />
  </figcaption>
</figure>

<!-- Dynamic axis scatter plot: defaults to field progress over time (CPS vs date
added) with a running-best line showing which releases moved the frontier -->
<h2>CPS Progress Over Time</h2>
<p>
  Each point is a model placed at its submission date; the dashed step line traces the
  running best ("SOTA frontier") CPS v1, so its jumps mark the models that set a new
  record when they were added. Use the axis/color/size selectors to compare models across
  any pair of metrics and parameters. The plot shows the same model cohort as the metrics
  table above, following the active task preset and table filters.
</p>
<DynamicScatter
  models={ACTIVE_MODELS}
  model_filter={in_cohort}
  x_key={METADATA_COLS.benchmark_added.key}
  y_key={ALL_METRICS.CPS.key}
  show_pareto_frontier
/>

<Readme>
  {#snippet title()}{/snippet}
  {#snippet model_count()}
    {ACTIVE_MODELS.filter(in_cohort).length}
  {/snippet}

  {#snippet best_report()}
    {#if best_entry}
      {@const { model: best_model, value } = best_entry}
      <span id="best-report">
        <a href="/models/{best_model.model_key}">{best_model.model_name}</a> leads the
        {col_preset} view with the best {@html primary_metric.label} of {format_num(
          value,
          primary_metric.format ?? `.3`,
        )}
        {primary_metric.unit ?? ``}.
      </span>
    {/if}
  {/snippet}
</Readme>
<!-- landing-only announcement; the shared MD note (warning + dataset + how to submit)
lives in MdNote and also renders on the /tasks/md page -->
<blockquote>
  🆕 <strong>New task — Molecular Dynamics.</strong> Matbench Discovery now scores how
  faithfully MLIPs reproduce finite-temperature observables of ab-initio MD (AIMD): radial
  distribution functions, vibrational density of states, pressure distributions, and
  single-point energy/force RMSEs. Explore the new metrics on the
  <a href="/tasks/md">MD leaderboard</a>.
</blockquote>
<MdNote />
<KappaNote warning={false} />

<h2>GitHub Activity</h2>
<p>
  Development activity and community engagement of MLIP GitHub repos. Points are sized by
  number of contributors and colored by number of commits over the last year.
</p>
<GitHubActivityScatter github_data={github_activity_data} />

<style>
  h1 {
    margin-block: -1.2em 1em;
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
  :is(.toggle-row, .downloads) {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
  }
  .toggle-row {
    gap: 8pt;
    font-size: smaller;
  }
  .downloads {
    gap: 1ex;
    margin-block: 1ex;
  }
  .downloads .download-btn {
    padding: 1pt 6pt;
    font: inherit;
  }
  .downloads a.download-btn {
    padding-inline-start: 2pt;
    &:not(:hover) {
      color: var(--text-color);
    }
  }
  .export-error {
    color: #ff6b6b;
    margin-block: 0.5em 1em;
    flex-basis: 100%;
    background-color: color-mix(in oklab, #ff6b6b 10%, transparent);
    padding: 1em;
    border-radius: 4px;
    border-inline-start: 4px solid #ff6b6b;
  }
  /* Caption Radar Container Styles */
  figcaption.caption-radar-container {
    display: grid;
    grid-template-columns: 1fr max-content;
    align-items: start;
    gap: 1em;
    font-size: 0.9em;
  }
</style>
