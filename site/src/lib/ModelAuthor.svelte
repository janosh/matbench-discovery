<script lang="ts">
  import { type Author, Icon } from '$lib'

  interface Props {
    author: Author
  }
  let { author }: Props = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
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
{/if}
