<script lang="ts">
  import { TableColumnToggleMenu } from '$lib'
  import { Tooltip } from 'svelte-zoo'
  import type { HeatmapColumn } from './types'

  // Props for this component
  interface Props {
    show_energy_only?: boolean
    show_noncompliant?: boolean
    columns?: HeatmapColumn[]
    on_filter_change?: (show_energy: boolean, show_noncomp: boolean) => void | undefined
  }

  // Extract props with defaults
  let {
    show_energy_only = $bindable(false),
    show_noncompliant = $bindable(false),
    columns = $bindable([]),
    on_filter_change = undefined,
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
      <svg><use href="#icon-info" /></svg>
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
      <svg><use href="#icon-info" /></svg>
    </Tooltip>
  </label>

  <TableColumnToggleMenu bind:columns bind:column_panel_open />
</div>

<style>
  .table-controls {
    display: flex;
    align-items: center;
    gap: 1rem;
    flex-wrap: wrap;
  }
  .filter-option {
    display: flex;
    align-items: center;
    gap: 0.3em;
  }
  /* Fix for sub and sup tags */
  :global(:is(sub, sup)) {
    font-size: 0.7em;
  }
</style>
