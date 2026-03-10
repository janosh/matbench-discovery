<script lang="ts">
  import { Icon } from 'matterviz'
  import { tooltip as tip } from 'svelte-multiselect/attachments'

  interface OptionInfo {
    value: string
    label: string
    tooltip?: string
    link?: string
  }
  let { selected = $bindable(``), options = [] }: {
    selected: string // currently selected value
    options: OptionInfo[] // options to display, either a record or an array of tuples
  } = $props()
  const target = { target: `_blank`, rel: `noopener noreferrer` }
</script>

<div class="selection-toggle">
  {#each options as { value, label, tooltip, link } (value)}
    <button
      class:active={selected === value}
      onclick={() => (selected = value)}
      {@attach tip({ content: tooltip })}
    >
      {@html label}
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
    padding: 4px 8px;
    border: 1px solid var(--border);
    background: transparent;
    cursor: pointer;
  }
  .selection-toggle button:hover {
    background: var(--nav-bg);
  }
  .selection-toggle button.active {
    background: var(--nav-bg);
  }
</style>
