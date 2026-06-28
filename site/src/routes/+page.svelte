<script lang="ts">
  import { replaceState } from '$app/navigation'
  import { page } from '$app/state'
  import { DATASETS, DISCOVERY_SETS, MetricsTable, SelectToggle } from '$lib'
  import {
    DynamicScatter,
    GitHubActivityScatter,
    RadarChart,
    SotaTimeline,
  } from '$lib/plot'
  import {
    ALL_METRICS,
    DISCOVERY_METRICS,
    discovery_set_toggle_options,
    GEO_OPT_SYMMETRY_METRICS,
    MD_METRICS,
    METADATA_COLS,
  } from '$lib/labels'
  import { make_combined_filter } from '$lib/metrics'
  import { find_best_model, MODELS } from '$lib/models.svelte'
  import {
    generate_csv,
    generate_excel,
    generate_png,
    generate_svg,
    handle_export,
  } from '$lib/table-export'
  import type { DiscoverySet, ModelData, SortDir } from '$lib/types'
  import Readme from '$root/readme.md'
  import MdNote from '$routes/tasks/md/md-note.md'
  import KappaNote from '$routes/tasks/phonons/kappa-note.md'
  import { format_num, Icon } from 'matterviz'
  import { onMount } from 'svelte'
  import { tooltip } from 'svelte-multiselect/attachments'
  import type { Snapshot } from './$types'
  import github_activity_data from './models/mlip-github-activity.json'

  const n_wbm_stable_uniq_protos = 32_942
  const n_wbm_uniq_protos = DATASETS.WBM.n_materials

  // Column presets focus the table on one task's metrics. Headline cols (CPS, F1, RMSD,
  // κ_SRME) stay visible across presets; supplementary discovery cols (TPR/TNR/RMSE) stay
  // hidden but remain toggleable via the per-column controls.
  const headline_metric_labels = new Set(
    [ALL_METRICS.CPS, DISCOVERY_METRICS.F1, ALL_METRICS.RMSD, ALL_METRICS.κ_SRME].map(
      (col) => col.label,
    ),
  )
  const supplementary_hidden = new Set([`TPR`, `TNR`, `RMSE`])
  const metadata_labels = new Set(Object.values(METADATA_COLS).map((col) => col.label))
  const col_presets = {
    Discovery: Object.values(DISCOVERY_METRICS),
    Phonons: [ALL_METRICS.κ_SRE],
    'Geo Opt': Object.values(GEO_OPT_SYMMETRY_METRICS),
    MD: Object.values(MD_METRICS),
  }
  type ColPreset = keyof typeof col_presets
  const col_preset_names = Object.keys(col_presets) as ColPreset[]
  const is_col_preset = (value: string | null): value is ColPreset =>
    col_preset_names.includes(value as ColPreset)
  const preset_default_sorts: Record<ColPreset, { column: string; dir: SortDir }> = {
    Discovery: { column: ALL_METRICS.CPS.key, dir: `desc` },
    Phonons: { column: ALL_METRICS.κ_SRME.key, dir: `asc` },
    'Geo Opt': { column: ALL_METRICS.RMSD.key, dir: `asc` },
    MD: { column: MD_METRICS.md_combined_score.key, dir: `desc` },
  }
  let show_non_compliant = $state(true)
  let show_energy_only = $state(false)
  let show_compliant = $state(true)
  const col_preset_options = col_preset_names.map((name) => ({
    value: name,
    label: name,
    tooltip: `Focus the table on ${name} metrics`,
  }))

  let show_heatmap = $state(true)
  let export_error: string | null = $state(null)

  let col_preset = $state<ColPreset>(`Discovery`)
  let preset_metric_labels = $derived(
    new Set([
      ...headline_metric_labels,
      ...col_presets[col_preset].map((col) => col.label),
    ]),
  )
  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
  let sort = $state<{ column: string; dir: SortDir }>({
    ...preset_default_sorts.Discovery,
  })
  let auto_sort_enabled = $state(true)
  let custom_col_config = $state(false)
  let previous_col_preset: ColPreset = `Discovery`
  let url_ready = $state(false)
  const sortable_header_selector = `thead th[role="button"]`

  $effect(() => {
    if (col_preset === previous_col_preset) return
    previous_col_preset = col_preset
    custom_col_config = false
    if (!auto_sort_enabled) return
    sort = { ...preset_default_sorts[col_preset] }
  })

  function handle_table_event(event: Event) {
    const target = event.target
    if (!(target instanceof Element)) return
    if (target.closest(sortable_header_selector)) auto_sort_enabled = false
    const reset_btn = target.closest(`button[aria-label="Reset all columns to defaults"]`)
    if (reset_btn?.closest(`details.column-toggles`)) custom_col_config = false
  }

  function handle_table_keydown(event: KeyboardEvent) {
    if ([`Enter`, ` `].includes(event.key)) handle_table_event(event)
  }

  function disable_col_preset_url(event: Event) {
    const target = event.target
    if (target instanceof HTMLInputElement && target.closest(`details.column-toggles`)) {
      custom_col_config = true
    }
  }

  onMount(() => {
    const params = page.url.searchParams
    const param_preset = params.get(`preset`)
    const next_preset = is_col_preset(param_preset) ? param_preset : `Discovery`
    const next_sort = { ...preset_default_sorts[next_preset] }
    const param_sort_col = params.get(`sort`)
    if (param_sort_col) next_sort.column = param_sort_col
    const param_sort_dir = params.get(`dir`)
    if (param_sort_dir === `asc` || param_sort_dir === `desc`)
      next_sort.dir = param_sort_dir
    const default_sort = preset_default_sorts[next_preset]
    auto_sort_enabled =
      next_sort.column === default_sort.column && next_sort.dir === default_sort.dir
    sort = next_sort

    const param_set = params.get(`set`) as DiscoverySet
    discovery_set = DISCOVERY_SETS.includes(param_set) ? param_set : `unique_prototypes`
    show_energy_only = params.get(`energy_only`) === `1`
    show_non_compliant = params.get(`non_compliant`) !== `0`
    show_compliant = params.get(`compliant`) !== `0`
    col_preset = next_preset
    previous_col_preset = next_preset
    url_ready = true
  })

  // Sync table state back to URL query params after the initial URL read.
  $effect(() => {
    if (!url_ready) return

    const new_params = new URLSearchParams()
    if (!custom_col_config) new_params.set(`preset`, col_preset)
    if (discovery_set !== `unique_prototypes`) new_params.set(`set`, discovery_set)
    const default_sort = preset_default_sorts[col_preset]
    if (sort.column !== default_sort.column) new_params.set(`sort`, sort.column)
    if (sort.dir !== default_sort.dir) new_params.set(`dir`, sort.dir)
    if (show_energy_only) new_params.set(`energy_only`, `1`)
    if (!show_non_compliant) new_params.set(`non_compliant`, `0`)
    if (!show_compliant) new_params.set(`compliant`, `0`)

    const new_url =
      new_params.size > 0 ? `${location.pathname}?${new_params}` : location.pathname
    if (new_url !== `${location.pathname}${location.search}`) {
      replaceState(new_url, page.state)
    }
  })

  // Export state object for handle_export
  let export_state = $derived({ export_error, show_non_compliant, discovery_set })

  let best_model = $derived(
    find_best_model(MODELS, { show_non_compliant, show_compliant, discovery_set }),
  )
  // Landing-page cohort: the metrics-table filters (energy/compliance) plus a base
  // predicate of "has discovery data for the selected set" (applied first internally)
  let in_cohort = $derived(
    make_combined_filter(
      (model: ModelData) => {
        const discovery = model.metrics?.discovery
        return (
          discovery !== null &&
          typeof discovery === `object` &&
          Boolean(discovery[discovery_set])
        )
      },
      show_energy_only,
      show_compliant,
      show_non_compliant,
    ),
  )

  export const snapshot: Snapshot = {
    capture: () => ({
      discovery_set,
      col_preset,
      custom_col_config,
      sort,
      auto_sort_enabled,
      show_non_compliant,
      show_energy_only,
      show_compliant,
      show_heatmap,
    }),
    restore: (values) => {
      custom_col_config = values.custom_col_config ?? custom_col_config
      auto_sort_enabled = values.auto_sort_enabled ?? auto_sort_enabled
      sort = values.sort ?? sort
      discovery_set = values.discovery_set ?? discovery_set
      col_preset = values.col_preset ?? col_preset
      previous_col_preset = col_preset
      show_non_compliant = values.show_non_compliant ?? show_non_compliant
      show_energy_only = values.show_energy_only ?? show_energy_only
      show_compliant = values.show_compliant ?? show_compliant
      show_heatmap = values.show_heatmap ?? show_heatmap
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
    <SelectToggle bind:selected={col_preset} options={col_preset_options} />
  </div>
  <!-- the test-set selector only affects discovery metrics, so only show it in the
  Discovery preset where those columns are visible -->
  {#if col_preset === `Discovery`}
    <div class="toggle-row compact">
      <span>Test set:</span>
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
    onchangecapture={disable_col_preset_url}
    onkeydowncapture={handle_table_keydown}
  >
    <MetricsTable
      col_filter={(col) =>
        metadata_labels.has(col.label)
          ? col.visible !== false
          : preset_metric_labels.has(col.label) && !supplementary_hidden.has(col.label)}
      {discovery_set}
      bind:sort
      bind:show_energy_only
      bind:show_non_compliant
      bind:show_compliant
      bind:show_heatmap
      show_energy_only_toggle
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
      style="color: var(--text-color)"
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
      estimators if an ensemble was used. DAF = Discovery Acceleration Factor measures how
      many more stable materials a model finds compared to random selection from the test
      set. The unique structure prototypes in the WBM test set
      {#if n_wbm_uniq_protos}
        have a
        <code>{format_num(n_wbm_stable_uniq_protos / n_wbm_uniq_protos, `.1%`)}</code>
        rate of stable crystals, meaning the max possible DAF is
        <code>
          ({format_num(n_wbm_stable_uniq_protos)} / {format_num(n_wbm_uniq_protos)})^−1 ≈
          {format_num(n_wbm_uniq_protos / n_wbm_stable_uniq_protos)}
        </code>.
      {:else}
        have an unknown rate of stable crystals (WBM n_materials unavailable).
      {/if}
    </div>
    <!-- CPS weight controls -->
    <RadarChart size={260} />
  </figcaption>
</figure>

<!-- Dynamic axis scatter plot -->
<p>Compare models across different metrics and parameters:</p>
<DynamicScatter models={MODELS} model_filter={in_cohort} />

<Readme>
  {#snippet title()}{/snippet}
  {#snippet model_count()}
    {MODELS.filter(in_cohort).length}
  {/snippet}

  {#snippet best_report()}
    {#if best_model}
      {@const { model_name, model_key, repo, paper, metrics = {} } = best_model}
      {@const discovery_metrics =
        typeof metrics?.discovery === `object` ? metrics.discovery : null}
      {@const { F1, R2, DAF } = discovery_metrics?.[discovery_set] ?? {}}
      <span id="best-report">
        <a href="/models/{model_key}">{model_name}</a> (<a href={paper}>paper</a>,
        <a href={repo}>code</a>) achieves the highest F1 score of {F1}, R<sup>2</sup> of {R2}
        and a discovery acceleration factor (DAF) of {DAF}
        (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
        dummy discovery in the already enriched test set containing 16% stable materials).
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

<h2>Progress Over Time</h2>
<p>
  Benchmark progress since launch: each point is a model placed at its submission date.
  The dashed step line traces the running best ("SOTA frontier") of the selected metric;
  labeled models are those that set a new record when they were added. The plot respects
  the compliance toggles of the metrics table above.
</p>
<SotaTimeline model_filter={in_cohort} style="height: 500px" />

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
    place-content: center;
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
  div.toggle-row {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
    gap: 8pt;
  }
  div.toggle-row > span {
    font-weight: 300;
    letter-spacing: 0.03em;
    opacity: 0.7;
  }
  /* slimmer secondary test-set selector (overrides SelectToggle's 4px block padding) */
  div.toggle-row.compact :global(.selection-toggle button) {
    padding-block: 1px;
  }
  div.downloads {
    display: flex;
    flex-wrap: wrap;
    gap: 1ex;
    justify-content: center;
    margin-block: 1ex;
    align-items: center;
  }
  div.downloads .download-btn {
    padding: 1pt 6pt;
    font: inherit;
  }
  div.export-error {
    color: #ff6b6b;
    margin-top: 0.5em;
    flex-basis: 100%;
    background-color: color-mix(in oklab, #ff6b6b 10%, transparent);
    padding: 1em;
    border-radius: 4px;
    border-left: 4px solid #ff6b6b;
    margin-bottom: 1em;
  }
  /* Caption Radar Container Styles */
  figcaption.caption-radar-container {
    display: grid;
    grid-template-columns: 1fr max-content;
    align-items: start;
    gap: 1em;
    font-size: 0.9em;
    background-color: transparent;
  }
  figure#metrics-table :global(:is(sub, sup)) {
    font-size: 0.7em;
  }
</style>
