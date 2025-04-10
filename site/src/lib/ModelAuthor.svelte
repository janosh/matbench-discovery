<script lang="ts">
  import type { Author } from '$lib'

  interface Props {
    author: Author
  }
  let { author }: Props = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
  <span title={affiliation}>{name}</span>
  {#if email}
    <a aria-label="Email" href="mailto:{email}">
      <svg><use href="#icon-mail"></use></svg>
    </a>
  {/if}
  {#if orcid}
    <a aria-label="ORCID" href={orcid}>
      <svg><use href="#icon-orcid"></use></svg>
    </a>
  {/if}
  {#if url}
    <a aria-label="Website" href={url}>
      <svg><use href="#icon-globe"></use></svg>
    </a>
  {/if}
  {#if github}
    <a aria-label="GitHub" href={github}>
      <svg><use href="#icon-github"></use></svg>
    </a>
  {/if}
{/if}

<style>
  span {
    font-size: smaller;
  }
  svg:has(use:only-child) {
    width: 16px;
    vertical-align: middle;
  }
</style>
