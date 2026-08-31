<script lang="ts">
  // Right-click menu for leaderboard rows: wraps a HeatmapTable and resolves the clicked
  // row to its model via the model-page link every row carries. Links, buttons and inputs
  // keep the browser's own menu (open in new tab, copy link).
  import { goto } from '$app/navigation'
  import { comparison } from '$lib/model-comparison.svelte'
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
    const href = target
      .closest(`tbody tr`)
      ?.querySelector(`a[href^="/models/"]`)
      ?.getAttribute(`href`)
    if (!href) return
    event.preventDefault()
    event.stopPropagation() // capture phase: preempts HeatmapTable's column menu
    model_key = href.slice(`/models/`.length)
    const point = { x: event.clientX, y: event.clientY }
    // Chromium fires contextmenu before pointerup; opening now would let that pointerup
    // light-dismiss the menu it just opened
    setTimeout(() => (at = point), 0)
  }
</script>

<div style="display: contents" oncontextmenucapture={open_menu}>
  {@render children()}
</div>
<ActionMenu {actions} bind:at trigger="none" aria-label="Model row actions" />
