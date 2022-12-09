<script lang="ts">
  const figs = import.meta.glob('$root/figures/*.{png,svg,pdf}', {
    eager: true,
    as: 'url',
  })
</script>

<h1>Figures</h1>

{#each Object.entries(figs) as [alt, src]}
  {@const filename = alt.split('/').at(-1)?.split('.')[0]}
  <h3>{filename?.replaceAll('-', ' ')}</h3>
  {#if src.endsWith('.pdf')}
    <embed {src} {alt} width="100%" height="500px" />
  {:else}
    <img {src} {alt} />
  {/if}
{/each}

<style>
  img {
    width: 100%;
  }
</style>
