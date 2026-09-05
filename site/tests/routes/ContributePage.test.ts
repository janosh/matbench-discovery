import Page from '$routes/contribute/+page.svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { mount } from '../index'

describe(`Contribute Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders markdown content with correct structure`, () => {
    expect(document.querySelector(`h1`)?.textContent).toBe(
      `How to submit new models to Matbench Discovery`,
    )

    const headings = [...document.querySelectorAll(`h2`)]
    expect(headings.length).toBeGreaterThanOrEqual(2)

    const expected_headings = [`Installation`, `Troubleshooting`]
    expected_headings.forEach((heading: string) => {
      expect(
        headings.some((h2) => h2.textContent?.includes(heading)),
        `Missing heading: ${heading}`,
      ).toBe(true)
    })

    const install_section = get_heading_section(document.body, `Installation`)
    expect(install_section?.textContent).toContain(`pip install -e ./matbench-discovery`)

    expect(document.querySelectorAll(`pre`).length).toBeGreaterThanOrEqual(2)
    const pypi_link = document.querySelector(
      `a[href*='pypi.org/project/matbench-discovery']`,
    )
    expect(pypi_link).not.toBeNull()
    expect(pypi_link?.getAttribute(`href`)).toMatch(
      /^https:\/\/pypi\.org\/project\/matbench-discovery\/?$/,
    )
    expect(pypi_link?.textContent).not.toBe(``)
  })

  it(`renders troubleshooting section with support links`, () => {
    const trouble_section = get_heading_section(document.body, `Troubleshooting`)

    const github_issue_link = trouble_section?.querySelector(
      `a[href*="github.com"][href*="issues"]`,
    )
    expect(github_issue_link).not.toBeNull()
    expect(github_issue_link?.getAttribute(`href`)).toMatch(
      /^https:\/\/github\.com\/.*\/issues/,
    )

    // the section has no community links (slack/forum/etc), so it must instruct
    // users to file an issue instead
    expect(github_issue_link?.textContent).toBe(`Open an issue on GitHub`)
    expect(trouble_section?.textContent).toContain(`happy to help!`)
  })
})

// content between the h2 whose text contains `heading_text` and the next h2
function get_heading_section(element: Element, heading_text: string): Element | null {
  const target_heading = [...element.querySelectorAll(`h2`)].find((h) =>
    h.textContent?.includes(heading_text),
  )

  if (!target_heading) return null

  const content = document.createElement(`div`)
  let current_node = target_heading.nextElementSibling

  while (current_node && current_node.tagName !== `H2`) {
    content.append(current_node.cloneNode(true))
    current_node = current_node.nextElementSibling
  }

  return content
}
