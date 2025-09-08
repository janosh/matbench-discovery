<script lang="ts">
  import { page } from '$app/state'
  import type { HTMLAttributes } from 'svelte/elements'
  import type { Snippet } from 'svelte'

  interface Props extends HTMLAttributes<HTMLNavElement> {
    routes: (string | [string, string])[]
    children?: Snippet
  }
  let { routes, children, ...rest }: Props = $props()

  let is_current = $derived((path: string) => {
    if (path === page.url.pathname) return `page`
    if (path !== `/` && page.url.pathname.includes(path)) return `page`
    return undefined
  })
</script>

<nav {...rest}>
  {#each routes as route, idx (route)}
    {@const [title, href] = Array.isArray(route) ? route : [route, route]}
    {#if idx > 0}<strong>&bull;</strong>{/if}
    <a {href} aria-current={is_current(href)} class="link">{title}</a>
  {/each}
  {#if children}<strong>&bull;</strong>{/if}
  {@render children?.()}
</nav>

<style>
  nav {
    display: flex;
    gap: 1em 1ex;
    place-content: center;
    place-items: center;
    margin: 2em auto;
    max-width: 45em;
    flex-wrap: wrap;
  }
  nav > a {
    padding: 0 5pt;
  }
</style>
