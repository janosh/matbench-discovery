<script lang="ts">
  import type { Author } from '$lib'
  import Icon from '@iconify/svelte'

  interface Props {
    author: Author
  }
  let { author }: Props = $props()
</script>

{#if author}
  {@const { name, email, orcid, affiliation, url, github } = author}
  <span title={affiliation}>{name}</span>
  {#if email}
    <a href="mailto:{email}"><Icon icon="ion:ios-mail" inline /></a>
  {/if}
  {#if orcid}
    <a href={orcid}><Icon icon="fa-brands:orcid" inline /></a>
  {/if}
  {#if url}
    <a href={url}><Icon icon="ion:ios-globe" inline /></a>
  {/if}
  {#if github}
    <a href={github}><Icon icon="octicon:mark-github" inline /></a>
  {/if}
{/if}

<style>
  span {
    font-size: smaller;
  }
</style>
