<script lang="ts">
  import { MODELS } from '$lib/models.svelte'
  import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
  import { format_num } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'

  let { model_key }: { model_key: string } = $props()

  // Includes energy-only models hidden by default and tracks live score weights.
  let ranks = $derived(model_metric_ranks(model_key, MODELS, RANKED_METRICS))

  // Mix theme-aware link blue (best) through muted purple to red (worst).
  const rank_color = (rank: number, n_models: number): string => {
    const frac = n_models > 1 ? (rank - 1) / (n_models - 1) : 0
    return `color-mix(in oklab, var(--link-color) ${Math.round(100 * (1 - frac))}%, hsl(0, 65%, 45%))`
  }

  const chip_title = ({ metric, rank, n_models, value }: (typeof ranks)[number]) => {
    const value_str = `${format_num(value, metric.format ?? `.3`)}${metric.unit ? ` ${metric.unit}` : ``}`
    // tooltip renders with allow_html, so use <br/> (not \r) for the line break
    return `Ranked ${rank} of ${n_models} models with a ${metric.label} of ${value_str}.<br/>${metric.description}`
  }
</script>

{#if ranks.length > 0}
  <section class="rank-card">
    <span class="rank-card-label">Leaderboard ranks</span>
    {#each ranks as rank_entry (rank_entry.metric.key)}
      {@const { metric, rank, n_models } = rank_entry}
      <a
        href={metric.rank_href}
        title={chip_title(rank_entry)}
        {@attach tooltip({ allow_html: true })}
      >
        <span class="metric-label">{@html metric.label}</span>
        <strong style:color={rank_color(rank, n_models)}>#{rank}</strong>
        <small>/{n_models}</small>
      </a>
    {/each}
  </section>
{/if}

<style>
  .rank-card {
    display: flex;
    flex-wrap: wrap;
    align-items: baseline;
    justify-content: center;
    gap: 3pt 1.4em;
    margin: 1em auto;
  }
  :is(.rank-card-label, .metric-label) {
    font-size: 0.9em;
  }
  :is(.rank-card-label, .metric-label, a small) {
    color: var(--text-secondary);
  }
  a {
    display: inline-flex;
    align-items: baseline;
    color: var(--text-color);
  }
  a strong {
    font-size: smaller;
    margin-left: 4pt;
  }
  a:hover .metric-label {
    text-decoration: underline;
  }
</style>
