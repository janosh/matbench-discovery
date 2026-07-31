<script lang="ts">
  import type { Author } from '$lib'
  import { get_org_logo } from '$lib/labels'
  import { Icon } from 'svelte-widgets'
  import { Contact, GitHub, Globe, Orcid } from 'svelte-widgets/icons'
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
  let { name, email, orcid, affiliation, url, github } = $derived(author)
  let org_logo = $derived(
    show_affiliation && affiliation ? get_org_logo(affiliation) : undefined,
  )
  let author_links = $derived([
    [email ? `mailto:${email}` : undefined, `Email`, Contact],
    [orcid, `Orcid`, Orcid],
    [url, `Website`, Globe],
    [github, `GitHub`, GitHub],
  ] as const)
</script>

<span {...rest}>
  <span title={affiliation}>{name}</span>
  {#if show_affiliation && affiliation}&ensp;<small>{affiliation}</small>{/if}
  {#if show_affiliation && org_logo}&nbsp;<Logo logo={org_logo} />{/if}
  {#each author_links as [href, label, icon] (label)}
    {#if href}
      <a {href} aria-label={label}><Icon {icon} /></a>
    {/if}
  {/each}
</span>
