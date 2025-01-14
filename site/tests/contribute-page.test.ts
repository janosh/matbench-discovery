import { beforeEach, describe, expect, it } from 'vitest'
import Page from '../src/routes/contribute/+page.svelte'

describe(`Contribute Page`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
    new Page({ target: document.body })
  })

  it(`renders markdown content`, () => {
    // Check for main headings
    const main_heading = document.body.querySelector(`h1`)
    expect(main_heading?.textContent).toBe(
      `How to submit new models to Matbench Discovery`,
    )

    // Check for installation section
    const install_heading = Array.from(document.body.querySelectorAll(`h2`)).find((h2) =>
      h2.textContent?.includes(`Installation`),
    )
    expect(install_heading).toBeDefined()

    // Check for code blocks
    const code_blocks = document.body.querySelectorAll(`pre`)
    expect(code_blocks.length).toBeGreaterThan(0)

    // Check for links
    const pypi_link = document.body.querySelector(
      `a[href*='pypi.org/project/matbench-discovery']`,
    )
    expect(pypi_link).toBeDefined()
  })

  it(`renders troubleshooting section`, () => {
    const trouble_heading = Array.from(document.body.querySelectorAll(`h2`)).find((h2) =>
      h2.textContent?.includes(`Troubleshooting`),
    )
    expect(trouble_heading).toBeDefined()

    const github_issue_link = Array.from(document.body.querySelectorAll(`a`)).find(
      (link) => link.href.includes(`github.com`) && link.href.includes(`issues`),
    )
    expect(github_issue_link).toBeDefined()
  })
})
