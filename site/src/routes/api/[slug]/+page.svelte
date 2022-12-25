<script lang="ts">
  import { page } from '$app/stores'
  import { error } from '@sveltejs/kit'

  $: slug = $page.url.pathname.split(`/`).at(-1)

  const routes = import.meta.glob(`../*.md`, { eager: true })

  $: component = routes[`../${slug}.md`]

  if (slug && !component) {
    throw error(404, `Page '${slug}' not found`)
  }
</script>

<svelte:component this={component?.default} />
