<script lang="ts">
  import spg_sankeys from '$figs/spg-sankeys.jsonl'
  import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
  import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
  import { GeoOptMetricsTable, ModelSelect, MODELS } from '$lib'
  import { order_models } from '$lib/fig-helpers'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import { make_table_filters } from '$lib/models.svelte'
  import { bind_url_params } from '$lib/url-state.svelte'
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
  const model_by_key = new Map(MODELS.map((model) => [model.model_key, model]))
  const model_key_by_label = new Map(
    MODELS.map((model) => [model.model_name, model.model_key]),
  )
  const plot_label_by_key = new Map([
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

  const date_added_ms = (key: string): number =>
    Date.parse(model_by_key.get(key)?.date_added ?? ``) || 0

  // newest models first; plot-only models without a date_added sort last
  const selectable_options = [...plot_label_by_key]
    .toSorted(([key_1], [key_2]) => date_added_ms(key_2) - date_added_ms(key_1))
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
  let filtered_sym_ops_sorted = $derived(
    sym_ops_sorted.filter(({ label }) => {
      const model_key = resolve_model_key(label)
      return model_key ? selected_model_key_set.has(model_key) : false
    }),
  )
  let filtered_spg_sankeys = $derived(
    spg_sankeys.models.filter(({ key }) => selected_model_key_set.has(key)),
  )

  const filters = make_table_filters()

  const read_url_params = (params: URLSearchParams) => {
    model_selection.read(params)
    filters.read(params)
  }
  bind_url_params(read_url_params, () => [
    model_selection.url_entry,
    ...filters.url_entries,
  ])

  const n_min_relaxed_structures =
    min(MODELS, ({ metrics }) =>
      typeof metrics?.geo_opt === `string`
        ? undefined
        : metrics?.geo_opt?.[`symprec=1e-2`]?.n_structures,
    ) ?? Infinity
</script>

<GeoOptReadme>
  {#snippet geo_opt_metrics_table()}
    <section class="full-bleed">
      <GeoOptMetricsTable {filters} />
    </section>
  {/snippet}
  {#snippet min_relaxed_structures()}
    <span>{format_num(n_min_relaxed_structures)}</span>
  {/snippet}
  {#snippet struct_rmsd_cdf_models()}
    <ScatterPlot
      series={struct_rmsd_sorted.map(({ label, auc, x, y }) => ({
        x,
        y,
        label: `${label} · AUC=${auc}`,
        markers: `line` as const,
      }))}
      x_axis={{ label: `RMSD (unitless)`, range: [0, 0.05] }}
      y_axis={{ label: `Cumulative`, format: `.0%`, range: [0, 1] }}
      style="height: 420px"
    />
  {/snippet}
  {#snippet sym_ops_diff_bar()}
    <div class="plot-controls bleed-1400">
      <ModelSelect
        options={selectable_options}
        bind:selected={model_selection.selected}
      />
    </div>

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
</GeoOptReadme>

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
