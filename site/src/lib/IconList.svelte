<script lang="ts">
  interface Props {
    icons: { id?: string; src?: string; name: string }[] | undefined
    style?: string
  }
  let { icons: logos = $bindable([]), style = `` }: Props = $props()
</script>

{#each logos ?? [] as logo (logo.id ?? logo.src)}
  {#if logo.id}
    <span title={logo.name} {style}>
      <svg style="margin: 0; height: 1em; width: 1em;">
        <use href="#{logo.id}"></use>
      </svg>
    </span>
  {:else if logo.src}
    <img
      src={logo.src}
      alt="{logo.name} logo"
      title={logo.name}
      style="margin: 0; height: 1em; {style}"
    />
  {/if}
{/each}

<style>
  svg,
  img {
    filter: grayscale(100%);
    height: 1em;
    width: auto;
    vertical-align: middle;
    margin: 0 !important;
  }
</style>
