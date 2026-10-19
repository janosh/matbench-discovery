import Page from '#routes/contribute/+page.svelte'
import { expect, it } from 'vite-plus/test'
import { doc_query, mount } from '../index.js'

it(`renders install and troubleshooting sections with package and support links`, () => {
  mount(Page, { target: document.body })
  expect(document.querySelector(`h1`)?.textContent).toBe(
    `How to submit new models to Matbench Discovery`,
  )

  const install_section = get_heading_section(document.body, `Installation`)
  expect(install_section.textContent).toContain(
    `activated Python 3.14+ virtual environment`,
  )
  expect(install_section.textContent).toContain(`uv pip install -e ./matbench-discovery`)
  expect(document.querySelectorAll(`pre`).length).toBeGreaterThanOrEqual(2)
  const pypi_link = doc_query(`a[href*='pypi.org/project/matbench-discovery']`)
  expect(pypi_link.getAttribute(`href`)).toMatch(
    /^https:\/\/pypi\.org\/project\/matbench-discovery\/?$/,
  )
  expect(pypi_link.textContent).not.toBe(``)

  // the section has no community links (slack/forum/etc), so it must instruct users to
  // file an issue instead
  const trouble_section = get_heading_section(document.body, `Troubleshooting`)
  const github_issue_link = doc_query(
    `a[href*="github.com"][href*="issues"]`,
    trouble_section,
  )
  expect(github_issue_link.getAttribute(`href`)).toMatch(
    /^https:\/\/github\.com\/.*\/issues/,
  )
  expect(github_issue_link.textContent).toBe(`Open an issue on GitHub`)
  expect(trouble_section.textContent).toContain(`happy to help!`)
})

// content between the h2 whose text contains `heading_text` and the next h2
function get_heading_section(element: Element, heading_text: string): Element {
  const target_heading = [...element.querySelectorAll(`h2`)].find((heading) =>
    heading.textContent?.includes(heading_text),
  )
  if (!target_heading) throw new Error(`Missing heading: ${heading_text}`)

  const content = document.createElement(`div`)
  let current_node = target_heading.nextElementSibling

  while (current_node && current_node.tagName !== `H2`) {
    content.append(current_node.cloneNode(true))
    current_node = current_node.nextElementSibling
  }

  return content
}
