<script lang="ts">
  import type { DiscoverySet, Metric, ModelData } from '$lib'
  import {
    MetricScatter,
    MetricsTable,
    model_is_compliant,
    MODELS,
    RadarChart,
    SelectToggle,
  } from '$lib'
  import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
  import { generate_png, generate_svg, handle_export } from '$lib/html-to-img'
  import { ALL_METRICS, DISCOVERY_SET_LABELS, METADATA_COLS } from '$lib/metrics'
  import { CPS_CONFIG } from '$lib/models.svelte'
  import Readme from '$root/readme.md'
  import KappaNote from '$routes/kappa-note.md'
  import { pretty_num } from 'elementari'

  let n_wbm_stable_uniq_protos = 32_942
  let n_wbm_uniq_protos = 215_488

  let show_non_compliant: boolean = $state(true)
  let show_energy_only: boolean = $state(false)
  let show_combined_controls: boolean = $state(true)
  let export_error: string | null = $state(null)

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    ...Object.fromEntries(
      [...ALL_METRICS, ...METADATA_COLS].map((col) => [col.label, true]),
    ),
    TPR: false,
    TNR: false,
    RMSE: false,
  })

  let best_model = $derived(
    MODELS.reduce((best: ModelData, md: ModelData) => {
      const best_F1 = best.metrics?.discovery?.full_test_set?.F1 ?? 0
      const md_F1 = md.metrics?.discovery?.full_test_set?.F1 ?? 0
      if (
        (!best_F1 || md_F1 > best_F1) &&
        (show_non_compliant || model_is_compliant(md))
      ) {
        return md
      }
      return best
    }, {} as ModelData),
  )

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)

  // Initialize with CPS as default metric
  let selected_metric:
    | keyof (typeof DEFAULT_CPS_CONFIG)[`parts`]
    | (typeof DEFAULT_CPS_CONFIG)[`key`] = $state(DEFAULT_CPS_CONFIG.key)
  let cps_scatter: Metric = {
    path: ``,
    label: `Combined Performance Score`,
    svg_label: DEFAULT_CPS_CONFIG.label,
    range: [0, 1],
    better: `higher`,
    description: DEFAULT_CPS_CONFIG.description,
  }
  let selected_scatter = $derived(
    { ...DEFAULT_CPS_CONFIG.parts, cps: cps_scatter }[selected_metric] as Metric,
  )
</script>

<h1 style="line-height: 0; margin-block: -1.2em 1em;">
  <img src="/favicon.svg" alt="Logo" width="60px" /><br />
  Matbench Discovery
</h1>

<figure style="margin-top: 3em;" id="metrics-table">
  <SelectToggle
    bind:selected={discovery_set}
    options={Object.entries(DISCOVERY_SET_LABELS).map(
      ([value, { title, tooltip, link }]) => ({ value, label: title, tooltip, link }),
    )}
  />

  <section class="table-wrapper">
    <div>
      <MetricsTable
        col_filter={(col) => visible_cols[col.label] ?? true}
        model_filter={() => true}
        {discovery_set}
        {show_combined_controls}
        {show_energy_only}
        show_noncompliant={show_non_compliant}
        config={CPS_CONFIG}
        style="width: 100%;"
      />
    </div>
  </section>

  <div class="downloads">
    Download table as
    {#each [[`SVG`, generate_svg], [`PNG`, generate_png]] as const as [label, generate_fn] (label)}
      <button
        class="download-btn"
        onclick={handle_export(generate_fn, label, export_error, {
          show_non_compliant,
          discovery_set,
        })}
      >
        {label}
      </button>
    {/each}
    <a
      href="/rss.xml"
      class="download-btn"
      title="Subscribe to be notified of new models"
    >
      <svg><use href="#icon-rss" /></svg>
      RSS
    </a>
    {#if export_error}
      <div class="export-error">
        {export_error}
      </div>
    {/if}
  </div>

  <!-- Radar Chart and Caption Container -->
  <figcaption class="caption-radar-container">
    <div
      style="flex: 1; background-color: var(--light-bg); padding: 0.2em 0.5em; border-radius: 4px;"
    >
      The <strong>CPS</strong> (Combined Performance Score) is a metric that weights
      discovery performance (F1), geometry optimization quality (RMSD), and thermal
      conductivity prediction accuracy (κ<sub>SRME</sub>). Use the radar chart to adjust
      the importance of each component.
      <br /><br />
      Training size is the number of materials used to train the model. For models trained
      on DFT relaxations, we show the number of distinct frames in parentheses). In cases where
      only the number of frames is known, we report the number of frames as the training set
      size. <code>(N=x)</code> in the Model Params column shows the number of estimators
      if an ensemble was used. DAF = Discovery Acceleration Factor measures how many more
      stable materials a model finds compared to random selection from the test set. The
      unique structure prototypes in the WBM test set have a
      <code>{pretty_num(n_wbm_stable_uniq_protos / n_wbm_uniq_protos, `.1%`)}</code>
      rate of stable crystals, meaning the max possible DAF is
      <code
        >({pretty_num(n_wbm_stable_uniq_protos)} / {pretty_num(n_wbm_uniq_protos)})^−1 ≈
        {pretty_num(n_wbm_uniq_protos / n_wbm_stable_uniq_protos)}</code
      >.
    </div>

    <!-- Radar chart with CPS weight controls -->
    <RadarChart size={260} />

    <!-- Model Size vs Performance plot -->
    <section style="width: 100%; margin-block: 1em;">
      <SelectToggle
        bind:selected={selected_metric}
        options={[
          {
            value: DEFAULT_CPS_CONFIG.key,
            label: DEFAULT_CPS_CONFIG.label,
            tooltip: DEFAULT_CPS_CONFIG.name,
          },
          ...Object.entries(DEFAULT_CPS_CONFIG.parts).map(
            ([value, { label, description }]) => ({
              value,
              label,
              tooltip: description,
            }),
          ),
        ]}
      />
      <h3 style="margin-block: 1em;">
        {@html selected_scatter.label} vs Model Size
        <small style="font-weight: lighter;">
          ({selected_scatter.better} = better, fewer model params = better)
        </small>
      </h3>
      <MetricScatter
        models={MODELS}
        config={selected_metric === DEFAULT_CPS_CONFIG.key ? CPS_CONFIG : undefined}
        metric={selected_metric !== DEFAULT_CPS_CONFIG.key ? selected_metric : undefined}
        y_label={selected_scatter.svg_label ?? selected_scatter.label}
        y_lim={selected_scatter.range}
        style="width: 100%; height: 300px;"
      />
    </section>

    <!-- Time-based scatter plot -->
    <section style="width: 100%;">
      <h3>
        {@html selected_scatter.label} over time
        <small style="font-weight: lighter;">
          ({selected_scatter.better} = better)
        </small>
      </h3>
      <MetricScatter
        models={MODELS}
        metric={selected_metric === DEFAULT_CPS_CONFIG.key ? undefined : selected_metric}
        config={selected_metric === DEFAULT_CPS_CONFIG.key ? CPS_CONFIG : undefined}
        y_label={selected_scatter.svg_label ?? selected_scatter.label}
        x_property="date_added"
        x_label="Date"
        range={selected_scatter.range}
        style="width: 100%; height: 300px;"
        date_range={[new Date(2024, 6, 1), null]}
      />
    </section>
  </figcaption>
</figure>

<Readme>
  {#snippet title()}{/snippet}
  {#snippet model_count()}
    {MODELS.filter((md) => show_non_compliant || model_is_compliant(md)).length}
  {/snippet}

  {#snippet best_report()}
    {#if best_model}
      {@const { model_name, model_key, repo, paper, metrics = {} } = best_model}
      {@const { F1, R2, DAF } = metrics?.discovery?.[discovery_set] ?? {}}
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
<KappaNote />
{#await import(`$routes/landing-page-figs.md`) then LandingPageFigs}
  <LandingPageFigs.default />
{/await}

<style>
  figure {
    margin: 0;
    display: grid;
    gap: 1ex;
  }
  /* Table wrapper for full-width placement */
  .table-wrapper {
    /* Use negative margin technique for full width */
    width: calc(100vw - 40px);
    margin-left: calc(-50vw + 50% + 20px);
    display: flex;
    justify-content: center;
  }
  figcaption {
    font-size: 0.9em;
    padding: 2pt 6pt;
    background-color: rgba(255, 255, 255, 0.06);
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
    background-color: rgba(255, 255, 255, 0.1);
    padding: 0 6pt;
    border-radius: 4pt;
    font: inherit;
  }
  div.export-error {
    color: #ff6b6b;
    margin-top: 0.5em;
    flex-basis: 100%;
  }

  /* Caption Radar Container Styles */
  figcaption.caption-radar-container {
    display: flex;
    flex-wrap: wrap;
    align-items: flex-start;
    gap: 1em;
    background-color: transparent;
  }
  figure#metrics-table :global(:is(sub, sup)) {
    font-size: 0.7em;
  }
</style>
