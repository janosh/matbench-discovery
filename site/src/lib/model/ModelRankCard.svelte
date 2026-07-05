<script lang="ts">
  import { MODELS } from '$lib/models.svelte'
  import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
  import type { ModelData } from '$lib/types'
  import { format_num } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'

  let { model_key, models = MODELS }: { model_key: string; models?: ModelData[] } =
    $props()

  // ranks span the full active roster (compliant + non-compliant + energy-only, so
  // the cohort can exceed the default leaderboard view, which hides energy-only
  // models); CPS/CMDS ranks track the session's current weight configs
  let ranks = $derived(model_metric_ranks(model_key, models, RANKED_METRICS))

  // green (best) -> red (worst) by rank position within the field
  const rank_color = (rank: number, n_models: number): string => {
    const frac = n_models > 1 ? (rank - 1) / (n_models - 1) : 0
    return `hsl(${Math.round(120 * (1 - frac))}, 60%, 42%)`
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
        class="rank-chip"
        href={metric.rank_href}
        style:--rank-color={rank_color(rank, n_models)}
        title={chip_title(rank_entry)}
        {@attach tooltip({ allow_html: true })}
      >
        <span class="metric-label">{@html metric.label}</span>
        <strong>#{rank}</strong>
        <small>/{n_models}</small>
      </a>
    {/each}
  </section>
{/if}

<style>
  .rank-card {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
    gap: 1ex;
    margin: 1em auto;
  }
  .rank-card-label {
    font-size: 0.9em;
    color: var(--text-secondary);
  }
  .rank-chip {
    display: inline-flex;
    align-items: baseline;
    gap: 3pt;
    padding: 1pt 7pt;
    border-radius: 4pt;
    background: color-mix(in srgb, var(--rank-color) 14%, transparent);
    border: 1px solid color-mix(in srgb, var(--rank-color) 45%, transparent);
    color: var(--text-color);
    text-decoration: none;
  }
  .rank-chip:hover {
    background: color-mix(in srgb, var(--rank-color) 26%, transparent);
  }
  .rank-chip .metric-label {
    font-size: 0.9em;
  }
  .rank-chip strong {
    color: var(--rank-color);
    filter: brightness(1.25);
  }
  .rank-chip small {
    color: var(--text-secondary);
  }
</style>
