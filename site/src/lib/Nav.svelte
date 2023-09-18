<script lang="ts">
  import { page } from '$app/stores'

  export let routes: string[]
  export let style: string | null = null

  $: is_current = (path: string) => {
    if (path === $page.url.pathname) return `page`
    if (path !== `/` && $page.url.pathname.includes(path)) return `page`
    return undefined
  }
</script>

<nav {style}>
  {#each routes as route, idx}
    {@const [title, href] = Array.isArray(route) ? route : [route, route]}
    {#if idx > 0}<strong>&bull;</strong>{/if}
    <a {href} aria-current={is_current(href)} class="link">{title}</a>
  {/each}
</nav>

<style>
  nav {
    display: flex;
    gap: 1em 1ex;
    place-content: center;
    margin: 1em auto 3em;
    max-width: 45em;
    flex-wrap: wrap;
    font-size: 1.1em;
  }
  nav > a {
    padding: 0 5pt;
  }
</style>
