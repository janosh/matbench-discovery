import { type Label, TableControls } from '$lib'
import { ALL_TRAINING_SETS, make_table_filters } from '$lib/models.svelte'
import { OPENNESS_OPTIONS } from '$lib/url-state.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`TableControls`, () => {
  const sample_columns: Label[] = [
    { key: `model`, label: `Model`, description: `Model name`, visible: true },
    { key: `f1`, label: `F1`, description: `F1 Score`, visible: true },
    { key: `daf`, label: `DAF`, description: `DAF Score`, visible: true },
    { key: `rmse`, label: `RMSE`, description: `RMSE`, visible: false },
  ]

  const summary_for = (text: string): HTMLElement => {
    const summary = [...document.querySelectorAll(`details.filter-menu summary`)].find(
      (el) => el.textContent?.includes(text),
    )
    if (!summary) throw new Error(`No filter summary found containing: ${text}`)
    return summary as HTMLElement
  }

  it(`renders filter dropdowns and heatmap toggle`, () => {
    mount(TableControls, { target: document.body })

    expect(summary_for(`Training data`)).toBeDefined()
    expect(summary_for(`Openness`)).toBeDefined()
    const labels = [...document.querySelectorAll(`label`)].map((label) =>
      label.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    expect(labels).toContain(`Heatmap`)
  })

  it(`training-data dropdown lists all datasets with require/exclude checkboxes`, async () => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    const dropdown = summary_for(`Training data`).closest(`details`)
    const rows = dropdown?.querySelectorAll(`.filter-row`) ?? []
    expect(rows).toHaveLength(ALL_TRAINING_SETS.length)

    // check `require` for the first dataset: require-mode filter becomes active,
    // its checkbox checks, and the summary shows a count badge
    const first_row = rows[0]
    const [require_box, exclude_box] =
      first_row.querySelectorAll<HTMLInputElement>(`input`)
    require_box.click()
    await tick()
    expect(filters.training[ALL_TRAINING_SETS[0]]).toBe(`require`)
    expect(require_box.checked).toBe(true)
    expect(summary_for(`Training data (1)`)).toBeDefined()

    // checking `exclude` on the same dataset flips the mode (mutually exclusive)
    exclude_box.click()
    await tick()
    expect(filters.training[ALL_TRAINING_SETS[0]]).toBe(`exclude`)
    expect(require_box.checked).toBe(false)
    expect(exclude_box.checked).toBe(true)

    // clear-filters button resets everything
    doc_query<HTMLButtonElement>(`button.clear-filters`).click()
    await tick()
    expect(filters.n_active).toBe(0)
  })

  it(`require means "must include, others allowed", not exclusivity`, () => {
    const filters = make_table_filters()
    filters.training = { MPtrj: `require` }
    // model trained on MPtrj plus extra datasets still matches
    const targets = `EFS_G`
    expect(filters.matches({ training_set: [`MPtrj`, `OMat24`, `sAlex`], targets })).toBe(
      true,
    )
    expect(filters.matches({ training_set: [`OMat24`], targets })).toBe(false)
  })

  it(`applies the built-in Compliant preset (old compliant cohort in one click)`, async () => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    const dropdown = summary_for(`Presets`).closest(`details`)
    const compliant_btn = [
      ...(dropdown?.querySelectorAll<HTMLButtonElement>(`button.preset`) ?? []),
    ].find((btn) => btn.textContent?.trim() === `Compliant`)
    if (!compliant_btn) throw new Error(`Compliant preset button not found`)
    compliant_btn.click()
    await tick()

    expect(filters.openness).toStrictEqual([`OSOD`])
    // mirrors old model_is_compliant: OSOD + only MP-anchored training data
    const targets = `EFS_G`
    const compliant = { training_set: [`MPtrj`, `MP 2022`], openness: `OSOD`, targets }
    const extra_data = { training_set: [`MPtrj`, `OMat24`], openness: `OSOD`, targets }
    const closed = { training_set: [`MPtrj`], openness: `CSOD`, targets }
    expect(filters.matches(compliant)).toBe(true)
    expect(filters.matches(extra_data)).toBe(false)
    expect(filters.matches(closed)).toBe(false)
  })

  it(`saves, applies and deletes user presets via localStorage`, async () => {
    localStorage.removeItem(`metrics-table-filter-presets`)
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    filters.set_training(`OMat24`, `exclude`)
    const dropdown = summary_for(`Presets`).closest(`details`)
    const input = dropdown?.querySelector<HTMLInputElement>(`form input`)
    if (!input) throw new Error(`preset name input not found`)
    input.value = `no-omat`
    input.dispatchEvent(new Event(`input`, { bubbles: true }))
    await tick()
    dropdown
      ?.querySelector(`form`)
      ?.dispatchEvent(new Event(`submit`, { bubbles: true, cancelable: true }))
    await tick()

    const stored = JSON.parse(
      localStorage.getItem(`metrics-table-filter-presets`) ?? `{}`,
    )
    expect(stored[`no-omat`]).toStrictEqual({
      training: { OMat24: `exclude` },
      openness: [`OSOD`, `OSCD`, `CSOD`, `CSCD`],
      targets: { F: `require` },
      fs_mode: `any`,
    })

    // clear filters, then re-apply via the saved preset button
    filters.clear()
    const preset_btn = [
      ...(dropdown?.querySelectorAll<HTMLButtonElement>(`button.preset`) ?? []),
    ].find((btn) => btn.textContent?.trim() === `no-omat`)
    preset_btn?.click()
    await tick()
    expect(filters.training).toStrictEqual({ OMat24: `exclude` })

    dropdown
      ?.querySelector<HTMLButtonElement>(`button[aria-label="Delete preset no-omat"]`)
      ?.click()
    await tick()
    expect(
      JSON.parse(localStorage.getItem(`metrics-table-filter-presets`) ?? `{}`),
    ).toStrictEqual({})
    // deletion is reactive: the preset button disappears from the dropdown
    expect(
      [...(dropdown?.querySelectorAll(`button.preset`) ?? [])].map((btn) =>
        btn.textContent?.trim(),
      ),
    ).not.toContain(`no-omat`)
  })

  it(`apply drops stale dataset keys and invalid modes from presets`, () => {
    const filters = make_table_filters()
    filters.apply({
      // deleted-dataset key and garbage mode could come from stale localStorage
      training: {
        MPtrj: `require`,
        'Renamed Dataset': `exclude`,
        OMat24: `bogus` as `exclude`,
      },
      openness: [`OSOD`, `bogus` as `OSCD`],
    })
    expect(filters.training).toStrictEqual({ MPtrj: `require` })
    expect(filters.openness).toStrictEqual([`OSOD`])
  })

  it(`openness dropdown toggles values but never hides the last one`, async () => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    const dropdown = summary_for(`Openness`).closest(`details`)
    const boxes = dropdown?.querySelectorAll<HTMLInputElement>(`input`) ?? []
    expect(boxes).toHaveLength(OPENNESS_OPTIONS.length)
    expect([...boxes].every((box) => box.checked)).toBe(true)

    boxes[1].click() // hide OSCD
    await tick()
    expect(filters.openness).toStrictEqual([`OSOD`, `CSOD`, `CSCD`])
    expect(summary_for(`Openness (3/4)`)).toBeDefined()
  })

  it(`opens and closes column visibility panel`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: sample_columns },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    const details = doc_query<HTMLDetailsElement>(`.column-toggles`)
    expect(toggle_btn.textContent?.trim()).toBe(`Columns`)

    // Should start closed
    expect(details.open).toBe(false)

    // Open panel
    toggle_btn.click()
    expect(details.open).toBe(true)

    // Verify column menu is visible
    expect(
      doc_query(`.column-menu`).querySelectorAll(`input[type="checkbox"]`),
    ).toHaveLength(sample_columns.length)

    // Close panel
    toggle_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles column visibility checkboxes`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: [...sample_columns] },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()

    const column_menu = doc_query(`.column-menu`)
    expect(column_menu.getAttribute(`role`)).toBe(`group`)

    const column_checkboxes =
      column_menu.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)
    expect(column_checkboxes).toHaveLength(sample_columns.length)
    const checkbox_labels = [...column_menu.querySelectorAll(`label`)].map((label) =>
      label.textContent?.trim(),
    )
    expect(checkbox_labels).toStrictEqual(sample_columns.map((column) => column.label))

    // Test toggling first checkbox (Model column, initially visible)
    const [first_checkbox] = column_checkboxes
    expect(first_checkbox.checked).toBe(true)

    first_checkbox.click()
    expect(first_checkbox.checked).toBe(false)

    first_checkbox.click()
    // Verify checkbox state changed
    expect(first_checkbox.checked).toBe(true)
  })
})
