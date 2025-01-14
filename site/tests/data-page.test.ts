import { beforeEach, describe, expect, it } from 'vitest'
import Page from '../src/routes/data/data-files-direct-download.md'

describe(`Data Page`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
    new Page({ target: document.body })
  })

  it(`renders data files list`, () => {
    const data_files_list = document.body.querySelector(`.data-files-list`)
    expect(data_files_list).toBeDefined()
    expect(data_files_list?.tagName.toLowerCase()).toBe(`ol`)
    expect(data_files_list?.childElementCount).toBeGreaterThan(10)
    expect(data_files_list?.textContent).toContain(`wbm-summary.csv.gz`)
  })

  it(`renders data files with valid URLs`, () => {
    const file_links = document.body.querySelectorAll(`.data-files-list a`)
    expect(file_links.length > 0).toBe(true)

    for (const link of file_links) {
      const url = link.getAttribute(`href`)
      expect(url).toBeTruthy()
      // URLs should be either Figshare or valid file paths
      expect(
        url?.includes(`figshare.com`) || url?.includes(`materialsproject`),
        `Invalid URL: ${url}`,
      ).toBeTruthy()
    }
  })
})
