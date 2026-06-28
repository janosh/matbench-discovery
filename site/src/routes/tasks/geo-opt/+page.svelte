<script lang="ts">
  import spg_sankeys from '$figs/spg-sankeys.jsonl'
  import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
  import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
  import { GeoOptMetricsTable, MODELS } from '$lib'
  import { order_models } from '$lib/fig-helpers'
  import { min } from 'd3-array'
  import { format_num } from 'matterviz'
  import { BarPlot, Sankey, sankey_from_links, ScatterPlot } from 'matterviz/plot'
  import GeoOptReadme from './geo-opt-readme.md'

  // payload models arrive pre-styled (colors + discovery-F1-desc leaderboard order) from the
  // figure_payload plugin; re-rank the two that want a different order (struct-rmsd by AUC
  // desc, sym-ops by symmetry-op-diff sigma asc). spg sankeys keep the leaderboard order.
  const struct_rmsd_sorted = order_models(struct_rmsd_cdf.models, (mdl) => -mdl.auc)
  const sym_ops_sorted = order_models(sym_ops_diff.models, (mdl) => mdl.sigma)

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
      <GeoOptMetricsTable show_non_compliant />
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
    <div class="sym-ops-list">
      {#each sym_ops_sorted as { label, sigma, x, y } (label)}
        <figure>
          <figcaption>{label} (σ={sigma})</figcaption>
          <BarPlot
            series={[{ x, y, label }]}
            y_axis={{ scale_type: `arcsinh` }}
            show_controls={false}
            style="height: 120px"
          />
        </figure>
      {/each}
    </div>
  {/snippet}
</GeoOptReadme>

<ul>
  {#each spg_sankeys.models as { key, label, labels, source, target, value } (key)}
    {@const data = sankey_from_links(source, target, value, labels)}
    <li>
      <h3>{label}</h3>
      <Sankey {data} show_controls={false} style="height: 300px; width: 100%" />
    </li>
  {/each}
</ul>

<style>
  ul {
    padding: 0;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
    gap: 3em 2em;
  }
  ul li {
    list-style: none;
  }
  ul h3 {
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
</style>
