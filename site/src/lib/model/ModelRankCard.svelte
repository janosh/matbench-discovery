<script lang="ts">
  import { MODELS } from '$lib/models.svelte'
  import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
  import type { ModelData } from '$lib/types'
  import { format_num } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'

  let { model_key, models = MODELS }: { model_key: string; models?: ModelData[] } =
    $props()

  // ranks span the full active roster (all models incl. energy-only, so
  // the cohort can exceed the default leaderboard view, which hides energy-only
  // models); CPS/CMDS ranks track the session's current weight configs
  let ranks = $derived(model_metric_ranks(model_key, models, RANKED_METRICS))

  // link-blue (best) -> red (worst) by rank position within the field, mixed in oklab
  // so mid-field ranks render as muted purple. Reusing --link-color (theme-aware)
  // keeps the best-rank blue consistent with the page's link styling
  const rank_color = (rank: number, n_models: number): string => {
    const frac = n_models > 1 ? (rank - 1) / (n_models - 1) : 0
    const blue_pct = Math.round(100 * (1 - frac))
    return `color-mix(in oklab, var(--link-color) ${blue_pct}%, hsl(0, 65%, 45%))`
  }

  const chip_title = ({ metric, rank, n_models, value }: (typeof ranks)[number]) => {
    const unit = metric.unit ? ` ${metric.unit}` : ``
    const value_str = `${format_num(value, metric.format ?? `.3`)}${unit}`
    // tooltip renders with allow_html, so use <br/> (not \r) for the line break
    return `Ranked ${rank} of ${n_models} models with a ${metric.label} of ${value_str}.<br/>${
      metric.description ?? ``
    }`
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
    gap: 4pt;
    color: var(--text-color);
  }
  a:hover .metric-label {
    text-decoration: underline;
  }
</style>
