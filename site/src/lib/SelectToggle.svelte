<script lang="ts">
  import { Tooltip } from 'svelte-zoo'

  interface OptionInfo {
    value: string
    label: string
    tooltip?: string
    link?: string
  }
  interface Props {
    selected: string // currently selected value
    options: OptionInfo[] // options to display, either a record or an array of tuples
  }
  let { selected = $bindable(``), options = [] }: Props = $props()
</script>

<div class="selection-toggle">
  {#each options as { value, label: option_label, tooltip, link } (value)}
    <Tooltip text={tooltip} tip_style="z-index: 2; font-size: 0.8em;">
      <button class:active={selected === value} onclick={() => (selected = value)}>
        {@html option_label}
        {#if link}
          <a
            href={link}
            target="_blank"
            rel="noopener noreferrer"
            aria-label="Info"
            style="line-height: 1;"
            onclick={(event) => event.stopPropagation()}
          >
            <svg><use href="#icon-info" /></svg>
          </a>
        {/if}
      </button>
    </Tooltip>
  {/each}
</div>

<style>
  .selection-toggle {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 5pt;
    margin-bottom: 5pt;
  }

  .selection-toggle button {
    padding: 4px 8px;
    border: 1px solid rgba(255, 255, 255, 0.05);
    background: transparent;
    cursor: pointer;
  }

  .selection-toggle button:hover {
    background: rgba(255, 255, 255, 0.05);
  }

  .selection-toggle button.active {
    background: rgba(255, 255, 255, 0.1);
  }
</style>
