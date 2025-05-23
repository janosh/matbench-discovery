<script lang="ts">
  import { Tooltip } from 'svelte-zoo'
  import { click_outside } from 'svelte-zoo/actions'
  import type { Label } from './types'

  interface Props {
    columns: Label[]
    column_panel_open?: boolean
  }

  let { columns = $bindable([]), column_panel_open = $bindable(false) }: Props = $props()

  function toggle_column_visibility(idx: number, event: Event) {
    columns[idx].visible = (event.target as HTMLInputElement).checked
    columns = [...columns] // needed to trigger reactivity
  }
</script>

<details
  class="column-toggles"
  bind:open={column_panel_open}
  use:click_outside={{ callback: () => (column_panel_open = false) }}
>
  <summary>
    Columns <svg><use href="#icon-columns"></use></svg>
  </summary>
  <div class="column-menu">
    {#each columns as col, idx (col.label + col.group + col.visible + idx)}
      <Tooltip style="white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
        {#snippet tip()}
          {#if col.description}
            {@html col.description}
          {/if}
        {/snippet}
        <label>
          <input
            type="checkbox"
            checked={col.visible !== false}
            onchange={(event) => toggle_column_visibility(idx, event)}
          />
          {@html col.label}
        </label>
      </Tooltip>
    {/each}
  </div>
</details>

<style>
  .column-toggles {
    position: relative;
  }
  .column-toggles summary {
    background: rgba(255, 255, 255, 0.1);
    padding: 2pt 6pt;
    border-radius: 4pt;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 4px;
  }
  .column-toggles summary:hover {
    background: rgba(255, 255, 255, 0.15);
  }
  .column-toggles summary::-webkit-details-marker {
    display: none;
  }
  .column-menu {
    font-size: 1.1em;
    position: absolute;
    right: 0;
    top: calc(100% + 4pt);
    background: #1c1c1c;
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 4pt;
    padding: 3pt 5pt;
    min-width: 150px;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(135px, 1fr));
    z-index: 1; /* needed to ensure column toggle menu is above HeatmapTable header row */
  }
  .column-menu label {
    display: inline-block;
    margin: 1px 2px;
    border-radius: 3px;
    line-height: 1.3em;
    height: 1.3em;
  }
  .column-menu label:hover {
    background: rgba(255, 255, 255, 0.1);
  }
  details :global(:is(sub, sup)) {
    transform: translate(-3pt, 6pt);
    font-size: 0.7em;
  }
</style>
