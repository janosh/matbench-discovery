<script lang="ts">
  import DataReadme from '$root/data/wbm/readme.md'
  import { onMount } from 'svelte'

  const figs = import.meta.glob(`$root/data/wbm/*.{png,svg,pdf}`, {
    eager: true,
    as: `url`,
  })

  onMount(() => {
    for (const img of document.querySelectorAll(`img`)) {
      const src = img.getAttribute(`src`)
      if (figs[`../data/wbm/${src}`]) {
        img.src = figs[`../data/wbm/${src}`]
      }
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
