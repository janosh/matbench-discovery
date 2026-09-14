import DATASETS from '$data/datasets.yml'
import { arr_to_str } from '$lib'
import type { Dataset } from '$lib/types'
import Page from '$routes/data/[slug]/+page.svelte'
import pkg from '$site/package.json'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`Dataset Detail Page`, () => {
  // Cover optional metadata, versioned datasets, and minimal entries.
  it.each([
    `MP 2022`,
    `NOMAD`,
    `ELEMENTA`,
    `MPtrj`,
    `SMAX`,
    `OMol25`,
    `OCx24`,
    `COSMOSDataset`,
    `MAD-1.5`,
    `WBM`,
  ])(`renders %s dataset correctly`, (dataset_key) => {
    const dataset = DATASETS[dataset_key]
    mount(Page, { target: document.body, props: { data: { dataset } } })

    expect(document.querySelector(`h1`)?.textContent).toBe(dataset.name)
    expect(document.querySelectorAll(`h1`)).toHaveLength(1)

    const meta_info = doc_query(`.meta-info`)
    expect(meta_info.textContent).toContain(`structures`)
    expect(meta_info.textContent?.toLowerCase()).toContain(`${dataset.access} access`)
    expect(meta_info.textContent?.toLowerCase()).not.toContain(dataset.role)
    expect(meta_info.textContent).toContain(dataset.license)
    for (const item of meta_info.children) {
      expect(item.querySelectorAll(`svg`), item.textContent).toHaveLength(1)
    }
    if (dataset.n_structures === null) {
      expect(meta_info.textContent).toContain(`Unknown number of structures`)
      expect(doc_query(`details`).textContent).toContain(
        `Experimental samples and electrodes are not atomic configurations.`,
      )
    }

    expect(doc_query(`.links`).querySelectorAll(`a`).length).toBeGreaterThan(0)
    expect(doc_query(`.description`).textContent).toMatch(/\S/)
    expect(document.querySelector(`.description p p`)).toBeNull()
    expect(
      document.querySelectorAll(
        `a[href="${pkg.repository}/blob/main/data/datasets.yml"]`,
      ),
    ).toHaveLength(2)
    if (dataset_key === `MPtrj`) {
      expect(doc_query(`#target-distributions`).textContent).toBe(`Target Distributions`)
      for (const count of [`1,580,395`, `7,944,833`, `49,295,660`, `14,223,555`]) {
        expect(document.body.textContent).toContain(count)
      }
    }
    if (dataset_key === `WBM`) {
      expect(document.body.textContent).toContain(`WBM Benchmark Details`)
      expect(document.body.textContent).toContain(`Downloading WBM data`)
      expect(document.body.textContent).toContain(`df_wbm.shape == (256_963, 18)`)
    }
    if (dataset_key === `ELEMENTA`) {
      const description = doc_query(`.description`)
      expect(description.textContent).toContain(
        `The full 210M-frame corpus is not public`,
      )
      expect(description.textContent).toContain(`coefficients sum to at most four`)
      expect(description.querySelector(`a`)?.getAttribute(`href`)).toBe(
        dataset.download_url,
      )
    }
  })

  // values used to be re-split on `:` after joining, truncating anything past a colon
  it.each([true, false])(`renders method parameters when present: %s`, (has_params) => {
    const dataset = DATASETS.WBM
    const params: Dataset[`params`] = {
      code: `VASP`,
      cutoff_energy: `520 eV: hard`,
      pseudopotentials: [`PBE`],
    }
    mount(Page, {
      target: document.body,
      props: {
        data: { dataset: { ...dataset, params: has_params ? params : undefined } },
      },
    })

    const items = [...doc_query(`.method-info`).querySelectorAll(`li`)].map((item) =>
      item.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    // scalar values go through arr_to_str (JSON-quoted), arrays are comma-joined
    expect(items).toEqual([
      `Method: ${arr_to_str(dataset.method)}`,
      ...(has_params
        ? [`Code: "VASP"`, `Cutoff Energy: "520 eV: hard"`, `Pseudopotentials: PBE`]
        : []),
    ])
  })
})
