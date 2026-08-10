import type { OrgLogo } from '$lib'
import { Magnetic, Meta } from 'svelte-widgets/icons'
import OrgLogos from '$lib/model/OrgLogos.svelte'
import { flushSync, type ComponentProps } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query, mount } from '../index'

describe(`OrgLogos.svelte`, () => {
  const mit_logo = {
    name: `Massachusetts Institute of Technology`,
    src: `/logos/massachusetts-institute-of-technology.svg`,
  } satisfies OrgLogo

  // mount, hover the preview, and return the resulting tooltip content element
  const open_tooltip = async (
    props: ComponentProps<typeof OrgLogos>,
  ): Promise<HTMLElement> => {
    mount(OrgLogos, { target: document.body, props })
    flushSync() // ensure the tooltip attachment effect has run before hovering
    doc_query(`.org-preview`).dispatchEvent(new PointerEvent(`pointerover`))
    return vi.waitFor(() => doc_query(`.custom-tooltip .tooltip-content`))
  }

  it(`renders one preview logo per org (src + icon types)`, () => {
    const org_logos: OrgLogo[] = [mit_logo, { name: `FAIR at Meta`, icon: Meta }]
    mount(OrgLogos, { target: document.body, props: { org_logos } })

    expect(document.querySelectorAll(`.org-preview .org-logo`)).toHaveLength(2)
    // src logo renders as <img>, icon logo as <span><svg>
    expect(doc_query<HTMLImageElement>(`img.org-logo`).src).toContain(mit_logo.src)
    expect(doc_query(`span.org-logo svg`).getAttribute(`aria-label`)).toBe(
      `FAIR at Meta logo`,
    )
    expect(document.querySelector(`.org-logo[title]`)).toBeNull()
  })

  it(`renders nothing when there are no org logos`, () => {
    mount(OrgLogos, { target: document.body, props: { org_logos: [] } })

    expect(document.querySelector(`.org-preview`)).toBeNull()
  })

  it(`shows grouped authors with escaped affiliations on hover`, async () => {
    const { innerHTML } = await open_tooltip({
      org_logos: [mit_logo, { name: `FAIR at Meta`, icon: Meta }],
      authors: [
        { name: `Ada Lovelace`, affiliation: `Massachusetts Institute of Technology` },
        { name: `Alan Turing`, affiliation: `Massachusetts Institute of Technology` },
        { name: `Grace Hopper`, affiliation: `University of Cambridge` },
        { name: `Yann LeCun`, affiliation: `FAIR at Meta` },
        { name: `Sam Houston`, affiliation: `Texas A&M University` },
      ],
    })

    // full affiliation names (not just logos) are shown
    expect(innerHTML).toContain(`Massachusetts Institute of Technology`)
    expect(innerHTML).toContain(`University of Cambridge`)
    // authors are grouped under their shared affiliation
    expect(innerHTML).toContain(`Ada Lovelace, Alan Turing`)
    expect(innerHTML).toContain(`Grace Hopper`)
    expect(innerHTML).toContain(`Yann LeCun`)
    expect(innerHTML).toContain(`Texas A&amp;M University`)
    // src affiliations render <img>, icon affiliations render inline <svg>
    expect(innerHTML).toContain(`<img`)
    expect(innerHTML).toContain(`<svg`)
  })

  it(`preserves stroked icon rendering in tooltip HTML`, async () => {
    const content_el = await open_tooltip({
      org_logos: [{ name: `Magnetic`, icon: Magnetic }],
    })
    const svg = doc_query(`svg`, content_el)

    expect(svg.getAttribute(`fill`)).toBe(`none`)
    expect(svg.getAttribute(`stroke`)).toBe(`currentColor`)
  })
})
