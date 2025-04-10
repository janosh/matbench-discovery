import Page from '$routes/data/data-files-direct-download.md'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'

describe(`Data Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders data files list with correct structure and content`, () => {
    const data_files_list = document.body.querySelector(`.data-files-list`)
    expect(data_files_list?.tagName.toLowerCase()).toBe(`ol`)

    // Verify proper number of data files are listed
    const list_items = data_files_list?.querySelectorAll(`li`)
    expect(list_items?.length).toBeGreaterThanOrEqual(16)

    // Check for expected data files in the list
    const list_text = Array.from(list_items || [])
      .map((item) => item.textContent || ``)
      .join(``)

    const essential_files = [
      `all_mp_tasks`,
      `mp_computed_structure_entries`,
      `mp_elemental_ref_entries`,
      `mp_energies`,
      `mp_patched_phase_diagram`,
      `mp_trj_json_gz`,
      `mp_trj_extxyz`,
      `wbm_computed_structure_entries`,
      `wbm_relaxed_atoms`,
      `wbm_initial_structures`,
      `wbm_initial_atoms`,
      `wbm_summary`,
      `wbm_dft_geo_opt_symprec_1e_2`,
      `wbm_dft_geo_opt_symprec_1e_5`,
      `phonondb_pbe_103_structures`,
      `phonondb_pbe_103_kappa_no_nac`,
    ]

    essential_files.forEach((filename) => {
      expect(list_text).toContain(filename)
    })

    // Verify each item has a link
    expect(list_items?.[0]?.querySelector(`a`)).toBeDefined()
  })

  it(`renders data files with valid and correctly formatted URLs`, () => {
    const file_links = document.body.querySelectorAll(`.data-files-list a`)
    expect(file_links.length).toBeGreaterThanOrEqual(16)

    // Verify all links have proper format
    for (const link of file_links) {
      const url = link.getAttribute(`href`)
      expect(url?.length).toBeGreaterThan(10)

      // Check URL has correct pattern
      const is_valid_url =
        /^https:\/\/(figshare\.com|.*materialsproject\.(org|com)|github\.com)\/.*$/.test(
          url || ``,
        )
      expect(is_valid_url, `Invalid URL format: ${url}`).toBe(true)

      // Link text should be meaningful
      expect(link.textContent?.trim().length).toBeGreaterThan(0)
    }
  })

  it(`renders file descriptions`, () => {
    // Extract descriptions from list items
    const descriptions = Array.from(document.body.querySelectorAll(`.data-files-list li`))
      .map((item) => {
        const link_text = item.querySelector(`a`)?.textContent || ``
        return item.textContent?.replace(link_text, ``).trim()
      })
      .filter((desc) => desc && desc.length > 0)

    // Verify descriptions exist and contain expected keywords
    expect(descriptions.length).toBeGreaterThan(5)
    expect(
      descriptions.some((desc) => /summary|structure|energy/i.test(desc || ``)),
    ).toBe(true)
  })
})
