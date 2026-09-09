<script lang="ts">
  import { ACTIVE_MODELS } from '$lib/models.svelte'
  import { model_metric_ranks, rank_color, RANKED_METRICS } from '$lib/rankings'
  import { format_num } from 'matterviz/labels'
  import { Popover } from 'svelte-widgets'

  let { model_key }: { model_key: string } = $props()

  // Rank against the active leaderboard cohort and track live score weights.
  let ranks = $derived(model_metric_ranks(model_key, ACTIVE_MODELS, RANKED_METRICS))
</script>

{#if ranks.length}
  <section class="rank-card">
    <span class="rank-card-label">Leaderboard ranks</span>
    {#each ranks as rank_entry (rank_entry.metric.key)}
      {@const { metric, rank, n_models, value } = rank_entry}
      <Popover trigger_mode="hover" trap_focus={false} aria-label="Metric rank">
        {#snippet trigger(trigger_props)}
          <a href={metric.rank_href} {...trigger_props}>
            <span class="metric-label">{@html metric.label}</span>
            <strong style:color={rank_color(rank, n_models)}>#{rank}</strong>
            <small>/{n_models}</small>
          </a>
        {/snippet}
        Ranked {rank} of {n_models} models with a {@html metric.label} of {format_num(
          value,
          metric.format ?? `.3`,
        )}{@html metric.unit ? ` ${metric.unit}` : ``}.<br />
        {@html metric.description ?? ``}
      </Popover>
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
