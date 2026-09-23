<script lang="ts">
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined-scores.svelte'
  import { ALL_METRICS } from '$lib/labels'
  import type { Label } from '$lib/types'
  import { format_num } from 'matterviz/labels'
  import { tooltip } from 'svelte-widgets/attachments'

  type WeightsConfig = Record<string, Label & { weight: number }>
  let {
    config = CPS_CONFIG,
    default_config = DEFAULT_CPS_CONFIG,
    title_label = ALL_METRICS.CPS,
  }: {
    config?: WeightsConfig
    default_config?: WeightsConfig
    title_label?: Label
  } = $props()

  const entries = $derived(Object.entries(config))
  const is_default = $derived(
    entries.every(([key, metric]) => metric.weight === default_config[key].weight),
  )
  const percentage = (weight: number) => format_num(weight * 100, `.1~f`)

  function change_weight(key: string, input: HTMLInputElement) {
    const value = input.valueAsNumber
    if (Number.isFinite(value) && value >= 0 && value <= 100) {
      const others = entries.filter(([other_key]) => other_key !== key)
      const total = others.reduce((sum, [, metric]) => sum + metric.weight, 0)
      config[key].weight = value / 100
      // A sole active component has no remaining proportions; split the rest equally.
      for (const [, metric] of others)
        metric.weight =
          (1 - value / 100) * (total > 0 ? metric.weight / total : 1 / others.length)
    }
    // Restore blank/invalid edits, including when the reactive value did not change.
    input.value = percentage(config[key].weight)
  }

  function reset_weights() {
    for (const [key, metric] of entries) metric.weight = default_config[key].weight
  }
</script>

<div class="score-weights" role="group" aria-label="{title_label.key} weights">
  <div class="heading">
    <strong title={title_label.description} {@attach tooltip()}>
      {@html title_label.label} weights
    </strong>
    <button type="button" disabled={is_default} onclick={reset_weights}>Reset</button>
  </div>
  <div class="inputs">
    {#each entries as [key, metric] (key)}
      <label title={metric.description} {@attach tooltip()}>
        <span>{@html metric.label}</span>
        <span class="value">
          <input
            type="number"
            min="0"
            max="100"
            step="any"
            aria-label="{metric.label.replace(/<[^>]*>/g, ``)} weight (%)"
            value={percentage(metric.weight)}
            onchange={(event) => change_weight(key, event.currentTarget)}
          />
          <span aria-hidden="true">%</span>
        </span>
      </label>
    {/each}
  </div>
  <small>Other weights adjust proportionally. Total: 100%.</small>
</div>

<style>
  .score-weights {
    max-width: 38rem;
    margin: 1.2em auto;
    font-size: 0.9rem;
  }
  .heading {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 0.5em;
    button {
      background: none;
      color: var(--link-color);
      padding: 0.1em 0.3em;
      font: inherit;
      &:disabled {
        opacity: 0.4;
      }
    }
  }
  .inputs {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(7.5rem, 1fr));
    gap: 0.5em;
  }
  label {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 0.4em;
    border: 1px solid var(--border);
    border-radius: 5px;
    padding: 0.35em 0.55em;
    &:focus-within {
      outline: 2px solid var(--link-color);
      outline-offset: 1px;
    }
  }
  .value {
    display: flex;
    align-items: baseline;
    color: var(--text-muted);
  }
  input {
    width: 5ch;
    min-width: 0;
    padding: 0;
    border: 0;
    outline: none;
    background: transparent;
    color: var(--text-color);
    font: inherit;
    font-variant-numeric: tabular-nums;
    text-align: right;
    appearance: textfield;
    &::-webkit-inner-spin-button,
    &::-webkit-outer-spin-button {
      appearance: none;
      margin: 0;
    }
  }
  small {
    display: block;
    margin-top: 0.5em;
    color: var(--text-muted);
    font-size: 0.85em;
  }
</style>
