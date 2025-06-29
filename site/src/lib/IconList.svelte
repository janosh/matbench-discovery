<script lang="ts">
  interface Props {
    icons: { id?: string; src?: string; name: string }[] | undefined
    [key: string]: unknown
  }
  let { icons: logos = $bindable([]), ...rest }: Props = $props()
</script>

{#each logos ?? [] as logo (logo.id ?? logo.src)}
  {#if logo.id}
    <span title={logo.name} {...rest}>
      <svg style="margin: 0; height: 1em; width: 1em">
        <use href="#{logo.id}"></use>
      </svg>
    </span>
  {:else if logo.src}
    {@const style = `margin: 0; height: 1em; ${rest.style ?? ``}`}
    <img
      src={logo.src}
      alt="{logo.name} logo"
      title={logo.name}
      {...rest}
      {style}
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
