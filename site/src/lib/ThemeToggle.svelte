<script lang="ts">
  import { Icon } from '$lib'
  import { onMount } from 'svelte'

  let theme: `light` | `dark` = $state(`dark`)

  onMount(() => {
    const prefers_dark = globalThis.matchMedia(`(prefers-color-scheme: dark)`).matches
    theme = localStorage.theme || (prefers_dark ? `dark` : `light`)
    document.documentElement.style.colorScheme = theme
  })

  function toggle_theme(): void {
    theme = ({ light: `dark`, dark: `light` } as const)[theme]
    document.documentElement.style.colorScheme = theme
    localStorage.theme = theme
  }
</script>

<button onclick={toggle_theme} title="Switch to {theme} theme">
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
