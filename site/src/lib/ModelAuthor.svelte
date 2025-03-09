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
      <iconify-icon icon="ion:ios-mail" inline></iconify-icon>
    </a>
  {/if}
  {#if orcid}
    <a aria-label="ORCID" href={orcid}>
      <iconify-icon icon="fa-brands:orcid" inline></iconify-icon>
    </a>
  {/if}
  {#if url}
    <a aria-label="Website" href={url}>
      <iconify-icon icon="ion:ios-globe" inline></iconify-icon>
    </a>
  {/if}
  {#if github}
    <a aria-label="GitHub" href={github}>
      <iconify-icon icon="octicon:mark-github" inline></iconify-icon>
    </a>
  {/if}
{/if}

<style>
  span {
    font-size: smaller;
  }
</style>
