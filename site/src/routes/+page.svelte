<script lang="ts">
  import { browser } from '$app/environment'
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import {
    DATASETS,
    DISCOVERY_SETS,
    DynamicScatter,
    GitHubActivityScatter,
    Icon,
    MetricsTable,
    RadarChart,
    SelectToggle,
  } from '$lib'
  import { ALL_METRICS, DISCOVERY_SET_LABELS, METADATA_COLS } from '$lib/labels'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import {
    generate_csv,
    generate_excel,
    generate_png,
    generate_svg,
    handle_export,
  } from '$lib/table-export'
  import type { DiscoverySet, ModelData, SortDir } from '$lib/types'
  import Readme from '$root/readme.md'
  import KappaNote from '$routes/tasks/phonons/kappa-note.md'
  import { format_num } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'
  import type { Snapshot } from './$types'
  import github_activity_data from './models/mlip-github-activity.json'

  const n_wbm_stable_uniq_protos = 32_942
  const n_wbm_uniq_protos = DATASETS.WBM.n_materials

  let table = $state({
    show_non_compliant: true,
    show_energy_only: false,
    show_combined_controls: true,
    show_compliant: true,
    show_heatmap: true,
  })
  let export_error: string | null = $state(null)
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
  let sort = $state({ column: `CPS`, dir: `desc` as SortDir })
  let url_initialized = false

  // Sync table state with URL query params (read on mount, write on change)
  $effect(() => {
    if (!browser) return
    const params = page.url.searchParams

    if (!url_initialized) {
      url_initialized = true
      const param_set = params.get(`set`) as DiscoverySet
      if (DISCOVERY_SETS.includes(param_set)) discovery_set = param_set
      const param_sort = params.get(`sort`)
      if (param_sort) sort.column = param_sort
      const param_dir = params.get(`dir`)
      if (param_dir === `asc` || param_dir === `desc`) sort.dir = param_dir
      if (params.get(`energy_only`) === `1`) table.show_energy_only = true
      if (params.get(`non_compliant`) === `0`) table.show_non_compliant = false
      if (params.get(`compliant`) === `0`) table.show_compliant = false
      return
    }

    const new_params = new URLSearchParams()
    if (discovery_set !== `unique_prototypes`) new_params.set(`set`, discovery_set)
    if (sort.column !== `CPS`) new_params.set(`sort`, sort.column)
    if (sort.dir !== `desc`) new_params.set(`dir`, sort.dir)
    if (table.show_energy_only) new_params.set(`energy_only`, `1`)
    if (!table.show_non_compliant) new_params.set(`non_compliant`, `0`)
    if (!table.show_compliant) new_params.set(`compliant`, `0`)

    const new_url = new_params.size ? `?${new_params}` : page.url.pathname
    if (new_url !== `${page.url.pathname}${page.url.search}`) {
      goto(new_url, { replaceState: true, keepFocus: true, noScroll: true })
    }
  })

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
      ) return md
      return best
    }, {} as ModelData),
  )

  export const snapshot: Snapshot = {
    capture: () => ({ discovery_set, table, sort }),
    restore: (values) => {
      discovery_set = values.discovery_set ?? discovery_set
      table = values.table ?? table
      sort = values.sort ?? sort
    },
  }
</script>

<h1>
  <img src="/favicon.svg" alt="Matbench Discovery Logo" width="60px" />
  Matbench Discovery
</h1>

<figure style="margin-top: 3em" id="metrics-table">
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
      bind:sort
      bind:show_energy_only={table.show_energy_only}
      bind:show_non_compliant={table.show_non_compliant}
      bind:show_compliant={table.show_compliant}
      bind:show_heatmap={table.show_heatmap}
    />
  </section>

  {#if export_error}
    <div class="export-error">
      {export_error}
    </div>
  {/if}

  <div class="downloads">
    Download table as
    {#each [
        [`SVG`, generate_svg],
        [`PNG`, generate_png],
        [`CSV`, generate_csv],
        [`Excel`, generate_excel],
      ] as const as
      [label, generate_fn]
      (label)
    }
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
      The training set column shows the number of materials used to train the model. For
      models trained on DFT relaxations, we show the number of distinct frames in
      parentheses. In cases where only the number of frames is known, we report the number
      of frames as the training set size. <code>(N=x)</code> in the Model Params column
      shows the number of estimators if an ensemble was used. DAF = Discovery Acceleration
      Factor measures how many more stable materials a model finds compared to random
      selection from the test set. The unique structure prototypes in the WBM test set
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
<DynamicScatter
  models={MODELS}
  model_filter={(model) => table.show_non_compliant || model_is_compliant(model)}
/>

<Readme>
  {#snippet title()}{/snippet}
  {#snippet model_count()}
    {MODELS.filter((md) => table.show_non_compliant || model_is_compliant(md)).length}
  {/snippet}

  {#snippet best_report()}
    {#if best_model}
      {@const { model_name, model_key, repo, paper, metrics = {} } = best_model}
      {@const discovery_metrics = typeof metrics?.discovery === `object`
      ? metrics.discovery
      : null}
      {@const { F1, R2, DAF } = discovery_metrics?.[discovery_set] ?? {}}
      <span id="best-report">
        <a href="/models/{model_key}">{model_name}</a> (<a href={paper}>paper</a>,
        <a href={repo}>code</a>) achieves the highest F1 score of {F1}, R<sup>2</sup> of {
          R2
        } and a discovery acceleration factor (DAF) of {DAF}
        (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
        dummy discovery in the already enriched test set containing 16% stable materials).
      </span>
    {/if}
  {/snippet}
</Readme>
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
