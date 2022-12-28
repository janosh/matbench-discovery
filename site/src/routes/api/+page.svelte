<script lang="ts">
  import { onMount } from 'svelte'

  onMount(() => {
    for (const img of [
      ...document.querySelectorAll(
        `img[src='https://img.shields.io/badge/-source-cccccc?style=flat-square']`
      ),
    ] as HTMLAnchorElement[]) {
      img.src = `https://img.shields.io/badge/source-blue?style=flat`
    }
  })
</script>

{#each Object.values(import.meta.glob(`./*.md`, { eager: true })) as file}
  <svelte:component this={file?.default} />
{/each}

<style>
  :global(hr) {
    border: none;
    margin: 3em;
  }
  :global(code) {
    line-height: 1em;
    border-radius: 4pt;
  }
</style>
