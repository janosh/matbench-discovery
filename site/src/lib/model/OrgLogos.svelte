<script lang="ts">
  import type { Author, OrgLogo } from '$lib/types'
  import { get_org_logo } from '$lib/labels'
  import { Popover } from 'svelte-widgets'
  import Logo from '../Logo.svelte'

  let {
    org_logos = [],
    authors = [],
  }: {
    org_logos?: OrgLogo[]
    authors?: Author[]
  } = $props()

  // Group authors by affiliation; fall back to supplied logos.
  let entries = $derived.by(() => {
    const groups: { logo?: OrgLogo; label: string; names: string[] }[] = []
    for (const { name, affiliation } of authors) {
      const label = affiliation || `Affiliation n/a`
      let group = groups.find((grp) => grp.label === label)
      if (!group) {
        const logo = affiliation ? get_org_logo(affiliation) : undefined
        group = { logo, label, names: [] }
        groups.push(group)
      }
      if (name) group.names.push(name)
    }
    if (groups.length > 0) return groups
    return org_logos.map((logo) => ({ logo, label: logo.name, names: [] as string[] }))
  })
</script>

{#if org_logos.length > 0}
  <Popover
    trigger_mode="hover"
    trap_focus={false}
    placement="left"
    aria-label="Authors and affiliations"
  >
    {#snippet trigger(trigger_props)}
      <span class="org-preview" role="button" tabindex="0" {...trigger_props}>
        {#each org_logos as logo (logo.name)}
          <Logo {logo} show_title={false} />
        {/each}
      </span>
    {/snippet}
    <div style="display: flex; flex-direction: column; gap: 5px; text-align: left">
      {#each entries as { logo, label, names } (label)}
        <div>
          <div style="display: flex; align-items: center; gap: 6px; font-weight: 600">
            {#if logo}<Logo {logo} show_title={false} />{/if}
            <span>{label}</span>
          </div>
          {#if names.length > 0}
            <div style="opacity: 0.7; font-size: 0.9em">{names.join(`, `)}</div>
          {/if}
        </div>
      {/each}
    </div>
  </Popover>
{/if}

<style>
  .org-preview {
    display: inline-flex;
    align-items: center;
    gap: var(--org-logo-gap, 0.3em);
    justify-content: center;
    font-size: 1.2em;
  }
</style>
