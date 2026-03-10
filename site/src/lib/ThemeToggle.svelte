<script lang="ts">
  import { Icon } from 'matterviz'
  import { onMount } from 'svelte'

  let theme: `light` | `dark` = $state(`dark`)
  let next = $derived(({ light: `dark`, dark: `light` } as const)[theme])
  let title = $derived(`Switch to ${next} theme`)

  onMount(() => {
    const prefers_dark = globalThis.matchMedia(`(prefers-color-scheme: dark)`).matches
    theme = localStorage.theme || (prefers_dark ? `dark` : `light`)
    document.documentElement.style.colorScheme = theme
    document.documentElement.dataset.theme = theme
  })

  function toggle_theme(): void {
    theme = next
    document.documentElement.style.colorScheme = theme
    document.documentElement.dataset.theme = theme
    localStorage.theme = theme
  }
</script>

<button onclick={toggle_theme} {title} aria-label={title} type="button">
  <Icon
    icon={({ light: `Sun`, dark: `Moon` } as const)[theme]}
    style="transform: scale(1.2)"
  />
</button>

<style>
  button {
    border-radius: 50%;
    display: flex;
    place-items: center;
    place-content: center;
    width: 2em;
    height: 2em;
  }
</style>
