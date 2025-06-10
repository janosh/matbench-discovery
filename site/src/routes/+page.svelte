<script lang="ts">
  import type { DiscoverySet, ModelData } from '$lib'
  import { DynamicScatter, MetricsTable, RadarChart, SelectToggle } from '$lib'
  import { CPS_CONFIG } from '$lib/combined_perf_score.svelte'
  import { ALL_METRICS, DISCOVERY_SET_LABELS, METADATA_COLS } from '$lib/labels'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import {
    generate_csv,
    generate_excel,
    generate_png,
    generate_svg,
    handle_export,
  } from '$lib/table-export'
  import Readme from '$root/readme.md'
  import KappaNote from '$routes/tasks/phonons/kappa-note.md'
  import { format_num } from 'elementari'
  import type { Snapshot } from './$types'

  let n_wbm_stable_uniq_protos = 32_942
  let n_wbm_uniq_protos = 215_488

  let table = $state({
    show_non_compliant: true,
    show_energy_only: false,
    show_combined_controls: true,
    show_compliant: true,
    show_heatmap: true,
  })
  let export_error: string | null = $state(null)

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    ...Object.fromEntries(
      [...Object.values(ALL_METRICS), ...Object.values(METADATA_COLS)].map((col) => [
        col.label,
        col.visible !== false,
      ]),
    ),
    TPR: false,
    TNR: false,
    RMSE: false,
  })

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)

  // Export state object for handle_export
  let export_state = $derived({
    export_error,
    show_non_compliant: table.show_non_compliant,
    discovery_set,
  })

  let best_model = $derived(
    MODELS.reduce((best: ModelData, md: ModelData) => {
      const best_discovery = best.metrics?.discovery
      const md_discovery = md.metrics?.discovery

      const best_F1_raw =
        (typeof best_discovery === `object` && best_discovery?.full_test_set?.F1) ?? 0
      const md_F1_raw =
        (typeof md_discovery === `object` && md_discovery?.full_test_set?.F1) ?? 0

      // Ensure F1 values are numbers
      const best_F1 = typeof best_F1_raw === `number` ? best_F1_raw : 0
      const md_F1 = typeof md_F1_raw === `number` ? md_F1_raw : 0

      if (
        (!best_F1 || md_F1 > best_F1) &&
        (table.show_non_compliant || model_is_compliant(md))
      )
        return md
      return best
    }, {} as ModelData),
  )

  export const snapshot: Snapshot = {
    capture: () => ({ discovery_set, table }),
    restore: (values) => ({ discovery_set, table } = values),
  }
</script>

<h1 style="line-height: 0; margin-block: -1.2em 1em;">
  <img src="/favicon.svg" alt="Logo" width="60px" /><br />
  Matbench Discovery
</h1>

<figure style="margin-top: 3em;" id="metrics-table">
  <SelectToggle
    bind:selected={discovery_set}
    options={Object.entries(DISCOVERY_SET_LABELS).map(
      ([value, { label, description: tooltip, link }]) => ({
        value,
        label,
        tooltip,
        link,
      }),
    )}
  />

  <section class="full-bleed">
    <MetricsTable
      col_filter={(col) => visible_cols[col.label] ?? true}
      model_filter={() => true}
      {discovery_set}
      bind:show_energy_only={table.show_energy_only}
      bind:show_non_compliant={table.show_non_compliant}
      bind:show_compliant={table.show_compliant}
      bind:show_heatmap={table.show_heatmap}
      cps_config={CPS_CONFIG}
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
    &emsp;Subscribe to new models via
    <a
      href="/rss.xml"
      class="download-btn"
      title="Subscribe to be notified of new models"
      style="color: var(--text-color);"
    >
      <svg><use href="#icon-rss" /></svg> RSS
    </a>
  </div>

  <figcaption class="caption-radar-container">
    <div
      style="background-color: var(--light-bg); padding: 0.2em 0.5em; border-radius: 4px;"
    >
      The <strong>CPS</strong> (Combined Performance Score) is a metric that weights
      discovery performance (F1), geometry optimization quality (RMSD), and thermal
      conductivity prediction accuracy (κ<sub>SRME</sub>). Use the radar chart to adjust
      the importance of each component.
      <br /><br />
      The training set column shows the number of materials used to train the model. For models
      trained on DFT relaxations, we show the number of distinct frames in parentheses. In
      cases where only the number of frames is known, we report the number of frames as the
      training set size. <code>(N=x)</code> in the Model Params column shows the number of
      estimators if an ensemble was used. DAF = Discovery Acceleration Factor measures how
      many more stable materials a model finds compared to random selection from the test
      set. The unique structure prototypes in the WBM test set have a
      <code>{format_num(n_wbm_stable_uniq_protos / n_wbm_uniq_protos, `.1%`)}</code>
      rate of stable crystals, meaning the max possible DAF is
      <code>
        ({format_num(n_wbm_stable_uniq_protos)} / {format_num(n_wbm_uniq_protos)})^−1 ≈
        {format_num(n_wbm_uniq_protos / n_wbm_stable_uniq_protos)}
      </code>.
    </div>
    <!-- CPS weight controls -->
    <RadarChart size={260} />
  </figcaption>
</figure>

<!-- Dynamic axis scatter plot -->
<p>Compare models across different metrics and parameters:</p>
<DynamicScatter
  models={MODELS}
  model_filter={(model) => table.show_non_compliant || model_is_compliant(model)}
  {discovery_set}
/>

<Readme>
  {#snippet title()}{/snippet}
  {#snippet model_count()}
    {MODELS.filter((md) => table.show_non_compliant || model_is_compliant(md)).length}
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
<KappaNote warning={false} />

<style>
  figure {
    margin: 0;
    display: grid;
    gap: 1ex;
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
    transition: background-color 0.2s ease;
  }
  div.downloads .download-btn:hover {
    background-color: rgba(255, 255, 255, 0.2);
  }
  div.export-error {
    color: #ff6b6b;
    margin-top: 0.5em;
    flex-basis: 100%;
    background-color: rgba(255, 107, 107, 0.1);
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
