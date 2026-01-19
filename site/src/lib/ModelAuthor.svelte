<script lang="ts">
  import { type Author } from '$lib'
  import { Icon } from 'matterviz'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    author,
    ...rest
  }: HTMLAttributes<HTMLSpanElement> & {
    author: Author
  } = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
  <span {...rest}>
    <small title={affiliation}>{name}</small>
    {#if email}
      <a aria-label="Email" href="mailto:{email}">
        <Icon icon="Contact" />
      </a>
    {/if}
    {#if orcid}
      <a aria-label="Orcid" href={orcid}>
        <Icon icon="Orcid" />
      </a>
    {/if}
    {#if url}
      <a aria-label="Website" href={url}>
        <Icon icon="Globe" />
      </a>
    {/if}
    {#if github}
      <a aria-label="GitHub" href={github}>
        <Icon icon="GitHub" />
      </a>
    {/if}
  </span>
{/if}
