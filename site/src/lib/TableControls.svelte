<script lang="ts">
  import { Tooltip } from 'svelte-zoo'

  // Props for this component
  interface Props {
    show_energy_only?: boolean
    show_noncompliant?: boolean
    visible_cols?: Record<string, boolean>
    on_filter_change?: (show_energy: boolean, show_noncomp: boolean) => void | undefined
    on_col_change?: (column: string, visible: boolean) => void | undefined
  }

  // Extract props with defaults
  let {
    show_energy_only = $bindable(false),
    show_noncompliant = $bindable(false),
    visible_cols = $bindable({}),
    on_filter_change = undefined,
    on_col_change = undefined,
  }: Props = $props()

  // Column panel state
  let column_panel_open = $state(false)

  // Handle filter checkbox changes
  function handle_noncompliant_change(event: Event) {
    const target = event.target as HTMLInputElement
    const checked = target.checked

    // Update both local state and call callback
    show_noncompliant = checked
    on_filter_change?.(show_energy_only, checked)
  }

  function handle_energy_only_change(event: Event) {
    const target = event.target as HTMLInputElement
    const checked = target.checked

    // Update both local state and call callback
    show_energy_only = checked
    on_filter_change?.(checked, show_noncompliant)
  }

  // Handle column visibility changes
  function handle_column_change(column: string, visible: boolean) {
    // Update local state
    visible_cols[column] = visible

    // Call callback
    on_col_change?.(column, visible)
  }
</script>

<div class="table-controls">
  <label class="filter-option">
    <input
      type="checkbox"
      checked={show_noncompliant}
      onchange={handle_noncompliant_change}
    />
    <span>Non-compliant models</span>
    <Tooltip>
      <span class="info-icon-small">ⓘ</span>
      {#snippet tip()}
        <span>
          Models can be non-compliant for multiple reasons:<br />
          - closed source (model implementation and/or train/test code)<br />
          - closed weights<br />
          - trained on more than the permissible training set (<a
            href="https://docs.materialsproject.org/changes/database-versions#v2022.10.28"
            >MP v2022.10.28 release</a
          >)<br />
          We still show these models behind a toggle as we expect them<br /> to nonetheless
          provide helpful signals for developing future models.
        </span>
      {/snippet}
    </Tooltip>
  </label>

  <label class="filter-option">
    <input
      type="checkbox"
      checked={show_energy_only}
      onchange={handle_energy_only_change}
    />
    <span>Energy-only models</span>
    <Tooltip text="Include models that only predict energy (no forces or stress)">
      <span class="info-icon-small">ⓘ</span>
    </Tooltip>
  </label>

  <details class="column-toggles" bind:open={column_panel_open}>
    <summary>
      Columns <span class="column-icon">⊞</span>
    </summary>
    <div class="column-menu">
      {#each Object.keys(visible_cols) as col (col)}
        <label>
          <input
            type="checkbox"
            checked={visible_cols[col]}
            onchange={(e) => handle_column_change(col, e.currentTarget.checked)}
          />
          {@html col}
        </label>
      {/each}
    </div>
  </details>
</div>

<style>
  .table-controls {
    display: flex;
    align-items: center;
    gap: 1rem;
    flex-wrap: wrap;
  }

  .column-toggles {
    position: relative;
    display: inline-block;
  }

  .column-toggles summary {
    background: rgba(255, 255, 255, 0.1);
    padding: 2pt 6pt;
    border-radius: 4pt;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 4px;
    white-space: nowrap;
  }

  .column-toggles summary:hover {
    background: rgba(255, 255, 255, 0.15);
  }

  .column-toggles summary::-webkit-details-marker {
    display: none;
  }

  .column-icon {
    font-size: 0.9em;
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
    z-index: 10;
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

  .info-icon-small {
    font-size: 0.8em;
    margin-left: 0.2em;
    opacity: 0.7;
    cursor: help;
  }

  .filter-option {
    display: flex;
    align-items: center;
    font-size: 0.85em;
    gap: 0.3em;
    cursor: pointer;
    white-space: nowrap;
  }

  /* Fix for sub and sup tags */
  :global(sub),
  :global(sup) {
    font-size: 0.7em;
  }
</style>
