<script lang="ts">
  import { type Author } from '$lib'
  import { Icon } from 'matterviz'
  import type { SVGAttributes } from 'svelte/elements'

  let { author, ...rest }: SVGAttributes<SVGSVGElement> & {
    author: Author
  } = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
  <small title={affiliation}>{name}</small>
  {#if email}
    <a aria-label="Email" href="mailto:{email}">
      <Icon icon="Contact" {...rest} />
    </a>
  {/if}
  {#if orcid}
    <a aria-label="Orcid" href={orcid}>
      <Icon icon="Orcid" {...rest} />
    </a>
  {/if}
  {#if url}
    <a aria-label="Website" href={url}>
      <Icon icon="Globe" {...rest} />
    </a>
  {/if}
  {#if github}
    <a aria-label="GitHub" href={github}>
      <Icon icon="GitHub" {...rest} />
    </a>
  {/if}
{/if}
