<script lang="ts">
  // Right-click menu for leaderboard rows: wraps a HeatmapTable. Links, buttons and inputs
  // keep the browser's own menu (open in new tab, copy link).
  import { goto } from '$app/navigation'
  import { comparison, row_model_key } from '$lib/model-comparison.svelte'
  import { MODELS } from '$lib/models.svelte'
  import type { Snippet } from 'svelte'
  import { ActionMenu } from 'svelte-widgets'
  import type { CmdAction } from 'svelte-widgets'

  let { children }: { children: Snippet } = $props()

  let at = $state<{ x: number; y: number } | null>(null)
  let model_key = $state(``)
  let model = $derived(MODELS.find((md) => md.model_key === model_key))
  let selected = $derived(comparison.keys.has(model_key))

  let actions: CmdAction[] = $derived.by(() => {
    if (!model) return []
    const { model_name } = model
    const n_selected = comparison.keys.size
    return [
      {
        id: `toggle`,
        label: selected
          ? `Remove ${model_name} from comparison`
          : `Add ${model_name} to comparison`,
        action: () => comparison.toggle(model_key),
      },
      {
        id: `open`,
        label:
          selected && n_selected > 1
            ? `Compare ${n_selected} selected models`
            : `Compare ${model_name} with…`,
        action: () => comparison.open_with(model_key),
      },
      {
        id: `page`,
        label: `Open ${model_name} model page`,
        action: () => void goto(`/models/${model_key}`),
      },
    ]
  })

  function open_menu(event: MouseEvent) {
    const target = event.target instanceof Element ? event.target : null
    if (!target || target.closest(`a, button, input, select`)) return
    const key = row_model_key(target.closest(`tbody tr`))
    if (!key) return
    event.preventDefault()
    event.stopPropagation() // capture phase: preempts HeatmapTable's column menu
    model_key = key
    at = { x: event.clientX, y: event.clientY }
  }
</script>

<div style="display: contents" oncontextmenucapture={open_menu}>
  {@render children()}
</div>
<!-- press dismissal (pointerdown outside or Escape) instead of the popover's native
light-dismiss: on macOS/Linux contextmenu fires on mousedown, so the same gesture's mouseup
would otherwise close the menu it just opened. Sized like the table rather than body prose;
svelte-widgets' menu padding (3pt) and item padding (2pt 6pt) defaults are overridden. -->
<ActionMenu
  {actions}
  bind:at
  trigger="none"
  dismiss={{ dismiss_on: `press` }}
  aria-label="Model row actions"
  style="font-size: 12px; --action-menu-padding: 0; --action-menu-item-padding: 3pt 8pt"
/>
