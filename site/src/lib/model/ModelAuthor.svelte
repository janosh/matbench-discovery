<script lang="ts">
  import type { Author } from '$lib'
  import { get_org_logo } from '$lib/labels'
  import { Icon } from 'matterviz'
  import type { HTMLAttributes } from 'svelte/elements'
  import Logo from '../Logo.svelte'

  let {
    author,
    show_affiliation = false,
    ...rest
  }: HTMLAttributes<HTMLSpanElement> & {
    author: Author
    show_affiliation?: boolean
  } = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
  {@const org_logo =
    show_affiliation && affiliation ? get_org_logo(affiliation) : undefined}
  <span {...rest}>
    <small title={affiliation}>{name}</small>
    {#if show_affiliation && affiliation}&ensp;{affiliation}{/if}
    {#if show_affiliation && org_logo}&nbsp;<Logo logo={org_logo} />{/if}
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
