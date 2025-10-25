<script lang="ts">
  import { Icon, type Label, TableColumnToggleMenu } from '$lib'
  import { tooltip } from 'svelte-multiselect/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    show_energy_only = $bindable(false),
    columns = $bindable([]),
    show_heatmap = $bindable(true),
    show_compliant = $bindable(true),
    show_non_compliant = $bindable(true),
    show_selected_only = $bindable(false),
    selected_count = 0,
    on_filter_change = undefined,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    show_energy_only?: boolean
    columns?: Label[]
    show_heatmap?: boolean
    show_compliant?: boolean
    show_non_compliant?: boolean
    show_selected_only?: boolean
    selected_count?: number
    on_filter_change?: (
      show_energy: boolean,
      show_non_compliant: boolean,
    ) => void | undefined
  } = $props()
</script>

<div class="table-controls" {...rest}>
  {#if selected_count > 0}
    <label>
      <input
        type="checkbox"
        bind:checked={show_selected_only}
        aria-label="Toggle between showing only selected models and all models"
      />
      {show_selected_only ? `Show all` : `Show only ${selected_count} selected`}
    </label>
  {/if}

  <label class="legend-item" title="Toggle visibility of compliant models">
    <span class="color-swatch" style="background-color: var(--compliant-color)"></span>
    <input
      type="checkbox"
      bind:checked={show_compliant}
      onchange={(evt) => {
        if (!(evt.target as HTMLInputElement).checked && !show_non_compliant) {
          show_non_compliant = true // Prevent hiding both compliant and non-compliant models
        }
      }}
    />
    Compliant models
  </label>

  <label class="legend-item" title="Toggle visibility of non-compliant models">
    <input
      type="checkbox"
      bind:checked={show_non_compliant}
      onchange={(evt) => {
        if (!(evt.target as HTMLInputElement).checked && !show_compliant) {
          show_compliant = true // Prevent hiding both compliant and non-compliant models
        }
      }}
    />
    <span
      class="color-swatch"
      style="background-color: var(--non-compliant-color)"
    ></span>
    Non-compliant models
    <span
      {@attach tooltip({
        content: `
      Models can be non-compliant for multiple reasons:<br />
      - closed source (model implementation and/or train/test code)<br />
      - closed weights<br />
      - trained on more than the permissible training set (<a
        href="https://docs.materialsproject.org/changes/database-versions#v2022.10.28"
      >MP v2022.10.28 release</a>)<br />
      We still show these models behind a toggle as we expect them<br />
      to nonetheless provide helpful signals for developing future models.`,
      })}
    >
      <Icon icon="Info" />
    </span>
  </label>
  <label>
    <input
      type="checkbox"
      checked={show_energy_only}
      onchange={(event: Event) => {
        const target = event.target as HTMLInputElement
        // Update both local state and trigger callback (if passed)
        show_energy_only = target.checked
        on_filter_change?.(target.checked, false)
      }}
    />
    Energy-only models
    <span
      title="Include models that only predict energy (no forces or stress)"
      {@attach tooltip()}
    >
      <Icon icon="Info" />
    </span>
  </label>

  <label>
    <input
      type="checkbox"
      bind:checked={show_heatmap}
      aria-label="Toggle heatmap colors"
    />
    Heatmap
  </label>

  <TableColumnToggleMenu bind:columns />
</div>

<style>
  div.table-controls {
    display: inline-flex;
    flex-wrap: wrap;
    justify-content: end;
    gap: 4pt 12pt;
    align-items: center;
    font-size: 2cqw;
  }
  label.legend-item {
    display: flex;
    align-items: center;
    gap: 0.3em;
  }
  span.color-swatch {
    width: 3pt;
    height: 20pt;
    border-radius: 1pt;
  }
</style>
