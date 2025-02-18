<script lang="ts">
  import Icon from '@iconify/svelte'
  import { click_outside } from 'svelte-zoo/actions'

  interface Props {
    visible_cols: Record<string, boolean>
    column_panel_open?: boolean
  }

  let { visible_cols = $bindable({}), column_panel_open = $bindable(false) }: Props =
    $props()
</script>

<details
  class="column-toggles"
  bind:open={column_panel_open}
  use:click_outside={{ callback: () => (column_panel_open = false) }}
>
  <summary>
    Columns <Icon icon="octicon:columns-16" inline />
  </summary>
  <div class="column-menu">
    {#each Object.keys(visible_cols) as col}
      <label>
        <input type="checkbox" bind:checked={visible_cols[col]} />
        {@html col}
      </label>
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
    position: absolute;
    right: 0;
    top: calc(100% + 4pt);
    background: #1c1c1c;
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 4pt;
    padding: 3pt 5pt;
    min-width: 150px;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(120px, 1fr));
  }
  .column-menu label {
    display: inline-block;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
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
