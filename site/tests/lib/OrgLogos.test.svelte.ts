import type { Author } from '$lib'
import type { OrgLogo } from '$lib/labels'
import OrgLogos from '$lib/model/OrgLogos.svelte'
import { flushSync } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`OrgLogos.svelte`, () => {
  // mount, hover the preview, and return the resulting tooltip content element
  const open_tooltip = async (props: {
    org_logos?: OrgLogo[]
    authors?: Author[]
  }): Promise<HTMLElement> => {
    mount(OrgLogos, { target: document.body, props })
    flushSync() // ensure the tooltip attachment effect has run before hovering
    doc_query(`.org-preview`).dispatchEvent(new MouseEvent(`mouseenter`))
    await new Promise((resolve) => setTimeout(resolve, 200))
    return doc_query(`.custom-tooltip .tooltip-content`)
  }

  it(`renders one preview logo per org (src + icon types)`, () => {
    const org_logos: OrgLogo[] = [
      { name: `Massachusetts Institute of Technology`, src: `/logos/mit.svg` },
      { name: `FAIR at Meta`, validated_icon: `LogoMeta` },
    ]
    mount(OrgLogos, { target: document.body, props: { org_logos } })

    expect(document.querySelectorAll(`.org-preview .org-logo`)).toHaveLength(2)
    // src logo renders as <img>, icon logo as <span><svg>
    expect(doc_query<HTMLImageElement>(`img.org-logo`).src).toContain(`/logos/mit.svg`)
    expect(document.querySelector(`span.org-logo svg`)).not.toBeNull()
  })

  it(`suppresses per-logo native titles to avoid competing tooltips`, () => {
    const org_logos: OrgLogo[] = [{ name: `MIT`, src: `/logos/mit.svg` }]
    mount(OrgLogos, { target: document.body, props: { org_logos } })

    expect(doc_query(`.org-logo`).hasAttribute(`title`)).toBe(false)
  })

  it.each([1, 2, 4])(`renders $n_logos inline logos`, (n_logos) => {
    const org_logos: OrgLogo[] = Array.from({ length: n_logos }, (_, idx) => ({
      name: `Org ${idx}`,
      src: `/logos/mit.svg`,
    }))
    mount(OrgLogos, { target: document.body, props: { org_logos } })

    const org_preview = doc_query(`.org-preview`)
    expect(org_preview.querySelectorAll(`.org-logo`)).toHaveLength(n_logos)
  })

  it(`renders nothing when there are no org logos`, () => {
    mount(OrgLogos, { target: document.body, props: { org_logos: [] } })

    expect(document.querySelector(`.org-preview`)).toBeNull()
  })

  it(`shows full org names and authors grouped by affiliation on hover`, async () => {
    const { innerHTML } = await open_tooltip({
      org_logos: [
        { name: `Massachusetts Institute of Technology`, src: `/logos/mit.svg` },
        { name: `FAIR at Meta`, validated_icon: `LogoMeta` },
      ],
      authors: [
        { name: `Ada Lovelace`, affiliation: `Massachusetts Institute of Technology` },
        { name: `Alan Turing`, affiliation: `Massachusetts Institute of Technology` },
        { name: `Grace Hopper`, affiliation: `University of Cambridge` },
        { name: `Yann LeCun`, affiliation: `FAIR at Meta` },
      ],
    })

    // full affiliation names (not just logos) are shown
    expect(innerHTML).toContain(`Massachusetts Institute of Technology`)
    expect(innerHTML).toContain(`University of Cambridge`)
    // authors are grouped under their shared affiliation
    expect(innerHTML).toContain(`Ada Lovelace, Alan Turing`)
    expect(innerHTML).toContain(`Grace Hopper`)
    expect(innerHTML).toContain(`Yann LeCun`)
    // src affiliations render <img>, icon affiliations render inline <svg>
    expect(innerHTML).toContain(`<img`)
    expect(innerHTML).toContain(`<svg`)
  })

  it(`escapes HTML special characters in affiliation names`, async () => {
    const content_el = await open_tooltip({
      org_logos: [{ name: `Texas A&M University`, src: `/logos/texas-a&m.svg` }],
      authors: [{ name: `Sam Houston`, affiliation: `Texas A&M University` }],
    })

    // ampersand survives as text content (escaped in the raw HTML string)
    expect(content_el.textContent).toContain(`Texas A&M University`)
    expect(content_el.innerHTML).toContain(`Texas A&amp;M University`)
  })
})
