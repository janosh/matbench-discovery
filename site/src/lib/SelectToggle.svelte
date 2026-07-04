<script lang="ts">
  import { Icon } from 'matterviz'
  import { Spinner } from 'matterviz/feedback'
  import { tooltip as tip } from 'svelte-multiselect/attachments'

  interface OptionInfo {
    value: string
    label: string
    tooltip?: string
    link?: string
    loading?: boolean // show a small spinner in the button (e.g. while its tab loads)
  }
  let {
    selected = $bindable(``),
    options = [],
  }: {
    selected: string // Currently selected value
    options: OptionInfo[] // Options to display, either a record or an array of tuples
  } = $props()
  const target = { target: `_blank`, rel: `noopener noreferrer` }
</script>

<div class="selection-toggle">
  {#each options as { value, label, tooltip, link, loading } (value)}
    <button
      class:active={selected === value}
      aria-pressed={selected === value}
      onclick={() => (selected = value)}
      {@attach tip({ allow_html: true, content: tooltip })}
    >
      {@html label}
      {#if loading}
        <Spinner
          style="--spinner-size: 0.9em; --spinner-border-width: 2px; --spinner-margin: 0"
        />
      {/if}
      {#if link}
        <a href={link} onclick={(event) => event.stopPropagation()} {...target}>
          <Icon icon="Info" style="transform: scale(1.2) translateY(-1px)" />
        </a>
      {/if}
    </button>
  {/each}
</div>

<style>
  .selection-toggle {
    display: flex;
    flex-wrap: wrap;
    place-content: center;
    gap: 8pt;
  }
  .selection-toggle button {
    display: inline-flex;
    align-items: center;
    gap: 0.5ex;
    padding: 4px 8px;
    border: 1px solid var(--border);
    background: transparent;
    cursor: pointer;
  }
  .selection-toggle button:hover {
    background: var(--nav-bg);
  }
  .selection-toggle button.active {
    border-color: var(--link-color);
    background: color-mix(in oklab, var(--link-color) 8%, var(--nav-bg));
    color: var(--link-color);
  }
</style>
