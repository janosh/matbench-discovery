<script lang="ts">
  import spg_sankeys from '$figs/spg-sankeys.jsonl'
  import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
  import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
  import GeoOptMetricsTable from '$lib/table/GeoOptMetricsTable.svelte'
  import ModelSelect from '$lib/ModelSelect.svelte'
  import { ACTIVE_MODELS, make_table_filters } from '$lib/models.svelte'
  import { by_benchmark_added_desc } from '$lib'
  import { order_models } from '$lib/fig-helpers'
  import {
    ALL_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    scatter_options_by_key,
  } from '$lib/labels'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import { bind_url_params, UrlPlotState } from '$lib/url-state.svelte'
  import { min } from 'd3-array'
  import { format_num } from 'matterviz/labels'
  import { pick_contrast_color } from 'matterviz/colors'
  import { BarPlot, Sankey, sankey_from_links, ScatterPlot } from 'matterviz/plot'
  import GeoOptReadme from './geo-opt-readme.md'

  // payload models arrive pre-styled (colors + discovery-F1-desc leaderboard order) from the
  // json_payload plugin; re-rank the two that want a different order (struct-rmsd by AUC
  // desc, sym-ops by symmetry-op-diff sigma asc). spg sankeys keep the leaderboard order.
  const struct_rmsd_sorted = order_models(struct_rmsd_cdf.models, (mdl) => -mdl.auc)
  const sym_ops_sorted = order_models(sym_ops_diff.models, (mdl) => mdl.sigma)
  const default_n_models = 5
  const plot = new UrlPlotState(
    {
      x: ALL_METRICS.RMSD.key,
      y: GEO_OPT_SYMMETRY_METRICS[`symmetry_match_1e-2`].key,
      sort: {
        column: ALL_METRICS.RMSD.key,
        dir: `asc`,
      },
    },
    scatter_options_by_key,
  )
  const model_by_key = new Map(ACTIVE_MODELS.map((model) => [model.model_key, model]))
  const plot_label_by_key = new Map([
    ...struct_rmsd_cdf.models.map(({ model_key, label }) => [model_key, label] as const),
    ...sym_ops_diff.models.map(({ model_key, label }) => [model_key, label] as const),
    ...spg_sankeys.models.map(({ model_key, label }) => [model_key, label] as const),
  ])

  // plot-only models missing from ACTIVE_MODELS get a null date so they sort last
  const undated_model = { dates: { benchmark_added: null } }
  const model_or_undated = (key: string) => model_by_key.get(key) ?? undated_model

  // newest models first; plot-only models without a benchmark date sort last
  const selectable_options = [...plot_label_by_key]
    .toSorted(([key_1], [key_2]) =>
      by_benchmark_added_desc(model_or_undated(key_1), model_or_undated(key_2)),
    )
    .map(([key, label]) => {
      const model_color = model_by_key.get(key)?.color ?? `gray`
      const text_color = pick_contrast_color({ background: model_color })
      return {
        label,
        value: key,
        style: {
          selected: `background: ${model_color}; color: ${text_color};`,
          option: ``,
        },
      }
    })

  const default_selected_keys = selectable_options
    .slice(0, default_n_models)
    .map((option) => String(option.value))

  const model_selection = new UrlModelSelection(() => ({
    options: selectable_options,
    defaults: default_selected_keys,
  }))
  let selected_model_key_set = $derived(new Set(model_selection.values))
  let filtered_struct_rmsd_sorted = $derived(
    struct_rmsd_sorted.filter(({ model_key }) => selected_model_key_set.has(model_key)),
  )
  let filtered_sym_ops_sorted = $derived(
    sym_ops_sorted.filter(({ model_key }) => selected_model_key_set.has(model_key)),
  )
  let filtered_spg_sankeys = $derived(
    spg_sankeys.models.filter(({ model_key }) => selected_model_key_set.has(model_key)),
  )

  const filters = make_table_filters()

  const read_url_params = (params: URLSearchParams) => {
    model_selection.read(params)
    filters.read(params)
    plot.read(params)
  }
  bind_url_params(read_url_params, () => [
    model_selection.url_entry,
    ...filters.url_entries,
    ...plot.url_entries,
  ])

  const n_min_relaxed_structures =
    min(
      ACTIVE_MODELS,
      ({ metrics }) => metrics?.geo_opt?.[`symprec=1e-2`]?.n_structures,
    ) ?? Infinity
</script>

<GeoOptReadme>
  {#snippet geo_opt_metrics_table()}
    <section class="full-bleed">
      <GeoOptMetricsTable {filters} bind:sort={plot.sort} />
    </section>
  {/snippet}
  {#snippet min_relaxed_structures()}
    <span>{format_num(n_min_relaxed_structures)}</span>
  {/snippet}
  {#snippet model_comparison_scatter()}
    <h3 id="metric-comparison">
      {@html scatter_axis_label(plot.y)} vs {@html scatter_axis_label(plot.x)}
    </h3>
    <p>
      The default view compares structure-matching RMSD (lower is better) with the
      fraction of matching spacegroups at <code>symprec=1e-2</code> (higher is better). Marker
      size defaults to model parameters and color to training-set size.
    </p>
    <DynamicScatter
      models={ACTIVE_MODELS}
      model_filter={(model) => model.metrics?.geo_opt != null}
      bind:x_key={plot.x}
      bind:y_key={plot.y}
      color_key={METADATA_COLS.n_training_materials.key}
      show_pareto_frontier
      style="height: 800px"
    />
  {/snippet}
  {#snippet diagnostic_model_picker()}
    <div class="plot-controls bleed-1400">
      <ModelSelect options={selectable_options} bind:value={model_selection.selected} />
    </div>
  {/snippet}
  {#snippet struct_rmsd_cdf_models()}
    {#if filtered_struct_rmsd_sorted.length > 0}
      <div
        class="rmsd-cdf"
        role="group"
        aria-label="RMSD CDF models: {filtered_struct_rmsd_sorted
          .map(({ label }) => label)
          .join(`, `)}"
      >
        <ScatterPlot
          series={filtered_struct_rmsd_sorted.map(({ label, auc, x, y }) => ({
            x,
            y,
            label: `${label} · AUC=${auc}`,
            markers: `line` as const,
          }))}
          x_axis={{ label: `RMSD (unitless)`, range: [0, 0.05] }}
          y_axis={{ label: `Cumulative`, format: `.0%`, range: [0, 1] }}
          style="height: 420px"
        />
      </div>
    {:else}
      <p class="empty-note">No models selected. Pick models above to compare.</p>
    {/if}
  {/snippet}
  {#snippet sym_ops_diff_bar()}
    <div class="sym-ops-list bleed-1400">
      {#each filtered_sym_ops_sorted as { label, sigma, x, y } (label)}
        <figure>
          <figcaption>{label} (σ={sigma})</figcaption>
          <BarPlot
            series={[{ x, y, label }]}
            y_axis={{ scale_type: `arcsinh` }}
            show_controls={false}
            style="height: 120px"
          />
        </figure>
      {:else}
        <p class="empty-note">No models selected. Pick models above to compare.</p>
      {/each}
    </div>
  {/snippet}
  {#snippet spg_sankeys()}
    <ul class="spg-sankeys bleed-1400">
      {#each filtered_spg_sankeys as { model_key, label, labels, source, target, value } (model_key)}
        {@const n_labels = labels.length}
        {@const data = sankey_from_links(
          source,
          target.map((target_idx) => target_idx + n_labels),
          value,
          [
            ...labels.map((spg_label) => `DFT ${spg_label}`),
            ...labels.map((spg_label) => `Relaxed ${spg_label}`),
          ],
        )}
        <li>
          <h3>{label}</h3>
          <Sankey {data} show_controls={false} style="height: 300px; width: 100%" />
        </li>
      {:else}
        <li class="empty-note">No models selected. Pick models above to compare.</li>
      {/each}
    </ul>
  {/snippet}
</GeoOptReadme>

<style>
  .plot-controls {
    display: flex;
    justify-content: center;
    margin-block: 1.5em;
  }
  .empty-note {
    grid-column: 1 / -1;
    text-align: center;
    opacity: 0.7;
  }
  .sym-ops-list {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 32rem), 1fr));
    gap: 2em;
  }
  .spg-sankeys {
    padding: 0;
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 26rem), 1fr));
    gap: 3em 2em;
    list-style: none;
  }
  .spg-sankeys h3 {
    text-align: center;
    margin: 0 0 0.5em;
  }
  .sym-ops-list figure {
    margin: 0;
  }
  .sym-ops-list figcaption {
    text-align: center;
    font-size: 0.9em;
  }
  @media (width >= 900px) {
    .sym-ops-list {
      grid-template-columns: repeat(2, minmax(0, 1fr));
    }
  }
  @media (width >= 1200px) {
    .spg-sankeys {
      grid-template-columns: repeat(3, minmax(0, 1fr));
    }
  }
</style>
