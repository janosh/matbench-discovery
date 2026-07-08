<script lang="ts">
  import type { Author } from '$lib'
  import { get_org_logo, type OrgLogo } from '$lib/labels'
  import { escape_html, ICON_DATA } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'
  import Logo from '../Logo.svelte'

  // org_logos = deduped, logo-matched affiliations shown as a compact preview strip.
  // authors = full author list used to build the richer hover tooltip (full org
  // names + which authors are affiliated where) without widening the table column.
  let {
    org_logos = [],
    authors = [],
  }: {
    org_logos?: OrgLogo[]
    authors?: Author[]
  } = $props()

  // Render a logo as a standalone HTML string (the tooltip content is injected via
  // innerHTML into document.body, so component-scoped styles don't apply — inline only).
  const logo_html = (logo: OrgLogo): string => {
    const style = `height: 1.1em; width: auto; flex: 0 0 auto; vertical-align: middle`
    if (logo.validated_icon && logo.validated_icon in ICON_DATA) {
      const { path, viewBox } = ICON_DATA[logo.validated_icon]
      const inner = path.trim().startsWith(`<`) ? path : `<path d="${path}" />`
      return `<svg viewBox="${viewBox}" fill="currentColor" style="${style}">${inner}</svg>`
    }
    if (logo.src) {
      return `<img src="${escape_html(logo.src)}" alt="" style="${style}; filter: grayscale(100%)" />`
    }
    return ``
  }

  // One tooltip entry per affiliation: its logo, full name, and the authors there.
  // Authors are grouped by affiliation (preserving first-seen order); falls back to
  // the matched org logos when no author metadata is available.
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

  let tooltip_content = $derived.by(() => {
    const rows_html = entries
      .map(({ logo, label, names }) => {
        const head = `<div style="display: flex; align-items: center; gap: 6px; font-weight: 600">${
          logo ? logo_html(logo) : ``
        }<span>${escape_html(label)}</span></div>`
        const author_names =
          names.length > 0
            ? `<div style="opacity: 0.7; font-size: 0.9em">${escape_html(names.join(`, `))}</div>`
            : ``
        return head + author_names
      })
      .join(`\n`)
    return `<div style="display: flex; flex-direction: column; gap: 5px; text-align: left">${rows_html}</div>`
  })
</script>

{#if org_logos.length > 0}
  <span
    class="org-preview"
    {@attach tooltip({ allow_html: true, content: tooltip_content, placement: `left` })}
  >
    {#each org_logos as logo (logo.name)}
      <Logo {logo} show_title={false} />
    {/each}
  </span>
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
