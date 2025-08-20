import Page from '$routes/contribute/+page.svelte'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'

describe(`Contribute Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders markdown content with correct structure`, () => {
    // Check main heading
    expect(document.body.querySelector(`h1`)?.textContent).toBe(
      `How to submit new models to Matbench Discovery`,
    )

    // Check section headings
    const headings = Array.from(document.body.querySelectorAll(`h2`))
    expect(headings.length).toBeGreaterThanOrEqual(2)

    const expected_headings = [`Installation`, `Troubleshooting`]
    expected_headings.forEach((heading: string) => {
      expect(
        headings.some((h2) => h2.textContent?.includes(heading)),
        `Missing heading: ${heading}`,
      ).toBe(true)
    })

    // Check installation section
    const install_section = get_heading_section(document.body, `Installation`)
    expect(install_section?.textContent).toContain(
      `Clone the repo and install matbench-discovery`,
    )

    // Check code blocks and PyPI link
    expect(document.body.querySelectorAll(`pre`).length).toBeGreaterThanOrEqual(2)
    const pypi_link = document.body.querySelector(
      `a[href*='pypi.org/project/matbench-discovery']`,
    )
    expect(pypi_link).toBeDefined()
    expect(pypi_link?.getAttribute(`href`)).toMatch(
      /^https:\/\/pypi\.org\/project\/matbench-discovery\/?$/,
    )
    expect(pypi_link?.textContent).toBeTruthy()
  })

  it(`renders API examples and code snippets`, () => {
    const code_content = Array.from(document.body.querySelectorAll(`pre`))
      .map((block) => block.textContent || ``)
      .join(`\n`)

    expect(code_content).toMatch(/import|from|def|class|matbench|pip/i)
  })

  it(`renders troubleshooting section with support links`, () => {
    const trouble_section = get_heading_section(document.body, `Troubleshooting`)

    // Check for GitHub issues link
    const github_issue_link = trouble_section?.querySelector(
      `a[href*="github.com"][href*="issues"]`,
    )
    expect(github_issue_link).toBeDefined()
    expect(github_issue_link?.getAttribute(`href`)).toMatch(
      /^https:\/\/github\.com\/.*\/issues/,
    )

    // Check for helpful text
    expect(trouble_section?.textContent).toMatch(/problem|issue|error|help|support/i)

    // If no community links, should have instructions for filing issues
    const has_community_links = Array.from(
      trouble_section?.querySelectorAll(`a`) || [],
    ).some(
      (link) =>
        link.textContent?.includes(`@`) ||
        [`slack`, `forum`, `gitter`].some((term) => link.href?.includes(term)),
    )

    if (!has_community_links) {
      expect(trouble_section?.textContent).toMatch(/file|open|create|new issue/i)
    }
  })
})

// Helper function to get a section following a heading with specific text
function get_heading_section(element: Element, heading_text: string): Element | null {
  const target_heading = Array.from(element.querySelectorAll(`h2`)).find((h) =>
    h.textContent?.includes(heading_text)
  )

  if (!target_heading) return null

  // Get content between this heading and the next h2
  const content = document.createElement(`div`)
  let current_node = target_heading.nextElementSibling

  while (current_node && current_node.tagName !== `H2`) {
    content.appendChild(current_node.cloneNode(true))
    current_node = current_node.nextElementSibling
  }

  return content
}
