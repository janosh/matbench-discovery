<script lang="ts">
  import Readme from '$root/readme.md'
</script>

<Readme />

<h2>Figures</h2>

{#each Object.entries(import.meta.glob( `$root/figures/*.{png,svg}`, { eager: true, as: 'url' } )) as [alt, src]}
  <img {src} {alt} />
{/each}

<style>
  img {
    width: 100%;
  }
</style>
