<script lang="ts">
  import spg_sankeys from '$figs/spg-sankeys.jsonl'
  import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
  import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
  import { GeoOptMetricsTable, ModelSelect, ACTIVE_MODELS } from '$lib'
  import { order_models } from '$lib/fig-helpers'
  import {
    ALL_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    scatter_options_by_key,
  } from '$lib/labels'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import { make_table_filters } from '$lib/models.svelte'
  import { DynamicScatter } from '$lib/plot'
  import type { SortState } from '$lib/url-state.svelte'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
  } from '$lib/url-state.svelte'
  import { min } from 'd3-array'
  import { format_num, pick_contrast_color } from 'matterviz'
  import { BarPlot, Sankey, sankey_from_links, ScatterPlot } from 'matterviz/plot'
  import GeoOptReadme from './geo-opt-readme.md'

  // payload models arrive pre-styled (colors + discovery-F1-desc leaderboard order) from the
  // json_payload plugin; re-rank the two that want a different order (struct-rmsd by AUC
  // desc, sym-ops by symmetry-op-diff sigma asc). spg sankeys keep the leaderboard order.
  const struct_rmsd_sorted = order_models(struct_rmsd_cdf.models, (mdl) => -mdl.auc)
  const sym_ops_sorted = order_models(sym_ops_diff.models, (mdl) => mdl.sigma)
  const default_n_models = 5
  const default_scatter_x = ALL_METRICS.RMSD.key
  const default_scatter_y = GEO_OPT_SYMMETRY_METRICS[`symmetry_match_1e-2`].key
  const default_table_sort: SortState = {
    column: ALL_METRICS.RMSD.key,
    dir: `asc`,
  }
  const model_by_key = new Map(ACTIVE_MODELS.map((model) => [model.model_key, model]))
  const model_key_by_label = new Map(
    ACTIVE_MODELS.map((model) => [model.model_name, model.model_key]),
  )
  const plot_label_by_key = new Map([
    ...struct_rmsd_cdf.models.map(
      ({ label }) => [model_key_by_label.get(label) ?? label, label] as const,
    ),
    ...sym_ops_diff.models.map(
      ({ label }) => [model_key_by_label.get(label) ?? label, label] as const,
    ),
    ...spg_sankeys.models.map(({ key, label }) => [key, label] as const),
  ])
  const selectable_model_keys = new Set(plot_label_by_key.keys())

  const resolve_model_key = (key_or_label: string): string | undefined => {
    const model_key = model_key_by_label.get(key_or_label) ?? key_or_label
    return selectable_model_keys.has(model_key) ? model_key : undefined
  }

  const benchmark_added_ms = (key: string): number =>
    Date.parse(model_by_key.get(key)?.dates.benchmark_added ?? ``) || 0

  // newest models first; plot-only models without a benchmark date sort last
  const selectable_options = [...plot_label_by_key]
    .toSorted(([key_1], [key_2]) => benchmark_added_ms(key_2) - benchmark_added_ms(key_1))
    .map(([key, label]) => {
      const model_color = model_by_key.get(key)?.color ?? `gray`
      const text_color = pick_contrast_color({ bg_color: model_color })
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
    from_url: resolve_model_key,
  }))
  let selected_model_key_set = $derived(new Set(model_selection.values))
  const is_selected_label = (label: string): boolean => {
    const model_key = resolve_model_key(label)
    return model_key ? selected_model_key_set.has(model_key) : false
  }
  let filtered_struct_rmsd_sorted = $derived(
    struct_rmsd_sorted.filter(({ label }) => is_selected_label(label)),
  )
  let filtered_sym_ops_sorted = $derived(
    sym_ops_sorted.filter(({ label }) => is_selected_label(label)),
  )
  let filtered_spg_sankeys = $derived(
    spg_sankeys.models.filter(({ key }) => selected_model_key_set.has(key)),
  )

  const filters = make_table_filters()
  let table_sort = $state({ ...default_table_sort })
  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)

  const read_url_params = (params: URLSearchParams) => {
    model_selection.read(params)
    filters.read(params)
    table_sort = sort_from_query(params, default_table_sort)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
  }
  bind_url_params(read_url_params, () => [
    model_selection.url_entry,
    ...filters.url_entries,
    ...sort_url_entries(table_sort, default_table_sort),
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
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
      <GeoOptMetricsTable {filters} bind:sort={table_sort} />
    </section>
  {/snippet}
  {#snippet min_relaxed_structures()}
    <span>{format_num(n_min_relaxed_structures)}</span>
  {/snippet}
  {#snippet model_comparison_scatter()}
    <h3>
      {@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}
    </h3>
    <p>
      The default view compares structure-matching RMSD (lower is better) with the
      fraction of matching spacegroups at <code>symprec=1e-2</code> (higher is better). Marker
      size defaults to model parameters and color to training-set size.
    </p>
    <DynamicScatter
      models={ACTIVE_MODELS}
      model_filter={(model) => model.metrics?.geo_opt != null}
      bind:x_key={scatter_x}
      bind:y_key={scatter_y}
      color_key={METADATA_COLS.n_training_materials.key}
      show_pareto_frontier
      style="height: 800px"
    />
  {/snippet}
  {#snippet diagnostic_model_picker()}
    <div class="plot-controls bleed-1400">
      <ModelSelect
        options={selectable_options}
        bind:selected={model_selection.selected}
      />
    </div>
  {/snippet}
  {#snippet struct_rmsd_cdf_models()}
    {#if filtered_struct_rmsd_sorted.length > 0}
      <div
        class="rmsd-cdf"
        role="group"
        aria-label={`RMSD CDF models: ${filtered_struct_rmsd_sorted.map(({ label }) => label).join(`, `)}`}
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
      {#each filtered_spg_sankeys as { key, label, labels, source, target, value } (key)}
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
  }
  .spg-sankeys li {
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
