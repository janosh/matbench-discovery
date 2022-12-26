<script lang="ts">
  import DataReadme from '$root/data/wbm/readme.md'
  import { onMount } from 'svelte'

  const figs = import.meta.glob(`$root/data/wbm/*.{png,svg,pdf}`, {
    eager: true,
    as: `url`,
  })

  onMount(() => {
    for (const img of document.querySelectorAll(`img`)) {
      const src = `../data/wbm/${img.getAttribute(`src`)}`
      if (figs[src]) img.src = figs[src]
    }
  })
</script>

<main>
  <DataReadme />
</main>

<style>
  :global(img) {
    max-width: 100%;
    margin: 1em auto;
  }
</style>
