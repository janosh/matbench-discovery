import TableControls from '$lib/table/TableControls.svelte'
import { fit_toolbar_links } from '$lib/table/fit-toolbar-links'
import type { Column } from 'matterviz/table'
import { ACTIVE_MODELS, ALL_TRAINING_SETS, make_table_filters } from '$lib/models.svelte'
import { OPENNESS_OPTIONS, type Openness } from '$lib/url-state.svelte'
import { comparison } from '$lib/model-comparison.svelte'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import { describe, expect, it, onTestFinished, vi } from 'vite-plus/test'
import { doc_query, filter_menu_trigger, mount, open_filter_menu } from '../index'

describe(`TableControls`, () => {
  it.each([755, 756, 1000])(
    `fits optional links to a %spx toolbar and updates after resize or label changes`,
    async (initial_width) => {
      document.body.innerHTML = `
        <div class="control-buttons" style="display: flex; column-gap: 12px">
          <div class="active-filters">Filters</div>
          <a>Submit a model</a>
          <a data-toolbar-optional>RSS</a>
          <button data-toolbar-optional>How to read the table</button>
          <div class="table-controls" style="font-size: 12px">Compare</div>
          <div style="display: contents"><button>Export</button></div>
        </div>`
      const toolbar = doc_query(`.control-buttons`)
      const controls = doc_query(`.table-controls`)
      const links = [...toolbar.querySelectorAll<HTMLElement>(`[data-toolbar-optional]`)]
      let width = initial_width
      let required_width = 756
      Object.defineProperty(toolbar, `clientWidth`, { get: () => width })
      vi.spyOn(toolbar, `getBoundingClientRect`).mockImplementation(() => {
        expect(toolbar.style.inlineSize).toBe(`max-content`)
        expect(controls.style.flexWrap).toBe(`nowrap`)
        expect(doc_query(`.active-filters`).style.display).toBe(`none`)
        return new DOMRect(0, 0, required_width, 20)
      })
      let resized = () => {}
      const resize = { observe: vi.fn(), unobserve: vi.fn(), disconnect: vi.fn() }
      vi.spyOn(globalThis, `ResizeObserver`).mockImplementation(
        class {
          observe = resize.observe
          unobserve = resize.unobserve
          disconnect = resize.disconnect
          constructor(callback: ResizeObserverCallback) {
            resized = () => callback([], resize)
          }
        },
      )
      const styled_elements = [toolbar, controls, doc_query(`.active-filters`)]
      const original_styles = styled_elements.map((element) => element.style.cssText)
      const cleanup = fit_toolbar_links(controls)
      expect(links.map((link) => link.hidden)).toEqual([width < 756, width < 756])
      expect(resize.observe).toHaveBeenCalledWith(toolbar)
      for (width of [755, 1000, 756]) {
        resized()
        expect(links.map((link) => link.hidden)).toEqual([width < 756, width < 756])
        expect(styled_elements.map((element) => element.style.cssText)).toEqual(
          original_styles,
        )
      }
      required_width = 806
      controls.textContent = `Compare (12)`
      await vi.waitFor(() => expect(links.every((link) => link.hidden)).toBe(true))
      required_width = 756
      controls.textContent = `Compare`
      await vi.waitFor(() => expect(links.some((link) => link.hidden)).toBe(false))
      cleanup?.()
      expect(resize.disconnect).toHaveBeenCalledOnce()
    },
  )

  const sample_columns: Column[] = [
    { id: `model`, label: `Model`, description: `Model name`, visible: true },
    { id: `f1`, label: `F1`, description: `F1 Score`, visible: true },
    { id: `daf`, label: `DAF`, description: `DAF Score`, visible: true },
    { id: `rmse`, label: `RMSE`, description: `RMSE`, visible: false },
  ]

  const mount_with_filters = async (models?: ModelData[], columns?: Column[]) => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters, models, columns } })
    await tick()
    return filters
  }

  it(`keeps desktop and mobile filter trees mounted and wires the sheet`, async () => {
    // both trees stay in the DOM (CSS toggles visibility) so SSR/client markup match
    const filters = await mount_with_filters()
    expect(document.querySelectorAll(`button.filter-menu-trigger`)).toHaveLength(3)
    const sheet_trigger = document.querySelector(`button.filter-sheet-trigger`)
    expect(sheet_trigger).toBeInstanceOf(HTMLButtonElement)
    if (!(sheet_trigger instanceof HTMLButtonElement)) return
    sheet_trigger.click()
    await tick()
    const dialog = doc_query<HTMLDialogElement>(`dialog[aria-label="Model filters"]`)
    expect(dialog.open).toBe(true)
    doc_query<HTMLInputElement>(
      `dialog[aria-label="Model filters"] input[aria-label="exclude OMat24"]`,
    ).click()
    await tick()
    expect(filters.training.OMat24).toBe(`exclude`)

    doc_query<HTMLButtonElement>(`button[aria-label="Close model filters"]`).click()
    await tick()
    expect(document.querySelector(`dialog[aria-label="Model filters"]`)).toBeNull()
  })

  it(`training-data dropdown lists all datasets with require/exclude checkboxes`, async () => {
    const filters = await mount_with_filters()

    const dropdown = await open_filter_menu(`Training data`)
    const boxes = dropdown?.querySelectorAll<HTMLInputElement>(`input`) ?? []
    expect(boxes).toHaveLength(2 * ALL_TRAINING_SETS.length)
    const require_boxes = [...boxes].filter((box) =>
      box.getAttribute(`aria-label`)?.startsWith(`require `),
    )
    const dataset_for = (box: HTMLInputElement): string =>
      box.getAttribute(`aria-label`)?.slice(`require `.length) ?? ``
    const usage_counts = require_boxes.map((box) => {
      const dataset = dataset_for(box)
      return ACTIVE_MODELS.filter(
        (model) => filters.matches(model) && model.training_sets.includes(dataset),
      ).length
    })
    expect(usage_counts).toStrictEqual(
      usage_counts.toSorted((count_left, count_right) => count_right - count_left),
    )

    // check `require` for the first dataset: require-mode filter becomes active,
    // its checkbox checks, and the summary shows a count badge
    const [require_box, exclude_box] = boxes
    const first_dataset = dataset_for(require_box)
    require_box.click()
    await tick()
    expect(filters.training[first_dataset]).toBe(`require`)
    expect(require_box.checked).toBe(true)
    expect(filter_menu_trigger(`Training data`).textContent).toContain(
      `Training data (1)`,
    )

    // checking `exclude` on the same dataset flips the mode (mutually exclusive)
    exclude_box.click()
    await tick()
    expect(filters.training[first_dataset]).toBe(`exclude`)
    expect(require_box.checked).toBe(false)
    expect(exclude_box.checked).toBe(true)

    // clear-filters button resets everything
    doc_query<HTMLButtonElement>(`button.clear-filters`).click()
    await tick()
    expect(filters.n_active).toBe(0)
  })

  it(`counts the current cohort with other constraints and exposes removable filter chips`, async () => {
    onTestFinished(() => comparison.keys.clear())
    const models = [
      [`facet-a`, [`MPtrj`, `OMat24`], `OSOD`, `EFS_G`],
      [`facet-b`, [`MPtrj`], `OSOD`, `EF_D`],
      [`facet-c`, [`OMat24`], `OSCD`, `EF_G`],
      [`facet-d`, [`MPtrj`, `OMat24`], `OSOD`, `E`],
    ].map(([model_key, training_sets, openness, targets]) => ({
      ...ACTIVE_MODELS[0],
      model_key,
      training_sets,
      openness,
      targets,
    })) as ModelData[]
    const filters = await mount_with_filters(models)
    const panel_text = async (name: string) =>
      (await open_filter_menu(name)).textContent?.replaceAll(/\s+/g, ` `)
    const remove_filter = async (label: string) => {
      doc_query<HTMLButtonElement>(`button[aria-label="Remove ${label} filter"]`).click()
      await tick()
    }
    expect(await panel_text(`Training data`)).toContain(`MPtrj (2)`)
    expect(await panel_text(`Targets`)).toContain(`forces (F) (3)`)
    await remove_filter(`require forces`)
    expect(filters.targets).toEqual({})
    expect(document.querySelector(`.active-filters`)).toBeNull()
    expect(await panel_text(`Training data`)).toContain(`MPtrj (3)`)

    filters.set_target(`F`, `require`)
    filters.set_training(`MPtrj`, `require`)
    filters.set_training(`OMat24`, `exclude`)
    filters.openness = [`OSOD`]
    await tick()
    // OMat24 ignores its own exclusion but retains MPtrj, required forces and OSOD.
    expect(await panel_text(`Training data`)).toContain(`OMat24 (1)`)
    expect(await panel_text(`Training data`)).toContain(`MPtrj (1)`)
    expect(await panel_text(`Openness`)).toContain(`OSOD (1)`)
    expect(await panel_text(`Openness`)).toContain(`OSCD (0)`)
    expect(await panel_text(`Targets`)).toContain(`stress (S) (0)`)
    await remove_filter(`exclude OMat24`)
    expect(filters.training).toEqual({ MPtrj: `require` })
    expect(await panel_text(`Targets`)).toContain(`stress (S) (1)`)

    filters.fs_mode = `gradient`
    await tick()
    expect(await panel_text(`Training data`)).toContain(`MPtrj (1)`)
    // Mode counts ignore gradient itself, retaining the other constraints.
    expect(await panel_text(`Targets`)).toContain(`direct (1)`)
    expect(await panel_text(`Targets`)).toContain(`gradient (1)`)
    await remove_filter(`Forces/stress: gradient`)
    await remove_filter(`Openness: OSOD`)
    await remove_filter(`require MPtrj`)
    expect(filters.fs_mode).toBe(`any`)
    expect(filters.openness).toEqual(OPENNESS_OPTIONS)
    expect(filters.training).toEqual({})

    comparison.keys.add(`facet-c`)
    filters.show_selected_only = true
    await tick()
    expect(await panel_text(`Training data`)).toContain(`MPtrj (0)`)
    expect(await panel_text(`Training data`)).toContain(`OMat24 (1)`)
    await remove_filter(`Selected models only`)
    expect(filters.show_selected_only).toBe(false)
    filters.set_target(`F`, `exclude`)
    filters.set_training(`MPtrj`, `require`)
    await tick()
    filters.set_training(`MPtrj`, `exclude`)
    await tick()
    await remove_filter(`exclude forces`)
    await remove_filter(`exclude MPtrj`)
    expect(filters.targets).toEqual({})
    expect(filters.training).toEqual({})
    expect(document.querySelector(`.active-filters`)).toBeNull()
    comparison.keys.clear()
  })

  it.each([false, true])(
    `applies Compliant models from training filters (mobile=%s)`,
    async (mobile) => {
      const filters = await mount_with_filters()
      if (mobile) {
        doc_query<HTMLButtonElement>(`button.filter-sheet-trigger`).click()
        await tick()
      }
      const dropdown = mobile
        ? doc_query(`dialog[aria-label="Model filters"]`)
        : await open_filter_menu(`Training data`)
      const compliant_btn = [
        ...(dropdown?.querySelectorAll<HTMLButtonElement>(`button`) ?? []),
      ].find((btn) => btn.textContent?.trim() === `Compliant models`)
      if (!compliant_btn) throw new Error(`Compliant preset button not found`)
      compliant_btn.click()
      await tick()

      expect(filters.openness).toStrictEqual([`OSOD`])
      // the preset selects OSOD models trained only on MP-anchored data
      const targets = `EFS_G`
      const filter_model = (training_sets: string[], openness: Openness) => ({
        training_sets,
        openness,
        targets,
      })
      expect(filters.matches(filter_model([`MPtrj`, `MP 2022`], `OSOD`))).toBe(true)
      expect(filters.matches(filter_model([`MPtrj`, `OMat24`], `OSOD`))).toBe(false)
      expect(filters.matches(filter_model([`MPtrj`], `CSOD`))).toBe(false)
      // Derived filter entries must follow in-place updates and deletion immediately.
      filters.set_target(`F`, `exclude`)
      expect(filters.matches(filter_model([`MPtrj`], `OSOD`))).toBe(false)
      filters.set_target(`F`, `exclude`)
      expect(filters.matches(filter_model([`MPtrj`], `OSOD`))).toBe(true)
    },
  )

  it(`openness dropdown toggles values but never hides the last one`, async () => {
    const filters = await mount_with_filters()

    const dropdown = await open_filter_menu(`Openness`)
    const boxes = [...(dropdown?.querySelectorAll<HTMLInputElement>(`input`) ?? [])]
    expect(boxes).toHaveLength(OPENNESS_OPTIONS.length)
    expect(boxes.every((box) => box.checked)).toBe(true)

    // uncheck all but the last; UI badge tracks the active count
    for (const box of boxes.slice(0, -1)) {
      box.click()
      await tick()
    }
    expect(filters.openness).toStrictEqual([OPENNESS_OPTIONS.at(-1)])
    expect(filter_menu_trigger(`Openness`).textContent).toContain(`Openness (1/4)`)

    // clicking the sole remaining option must not empty the filter
    boxes.at(-1)?.click()
    await tick()
    expect(filters.openness).toStrictEqual([OPENNESS_OPTIONS.at(-1)])
    expect(boxes.at(-1)?.checked).toBe(true)
  })

  it(`opens, updates and closes the columns panel with its heatmap setting`, async () => {
    const filters = make_table_filters()
    mount(TableControls, {
      target: document.body,
      props: { columns: [...sample_columns], filters },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    const details = doc_query<HTMLDetailsElement>(`.column-toggles`)
    expect(toggle_btn.textContent?.trim()).toBe(`Columns`)
    expect(details.open).toBe(false)

    toggle_btn.click()
    expect(details.open).toBe(true)

    const column_menu = doc_query(`.column-menu`)
    expect(column_menu.getAttribute(`role`)).toBe(`group`)
    const heatmap = doc_query<HTMLInputElement>(`[aria-label="Toggle heatmap colors"]`)
    expect(column_menu.contains(heatmap)).toBe(true)
    expect(heatmap.checked).toBe(true)
    heatmap.click()
    await tick()
    expect(filters.show_heatmap).toBe(false)
    expect(details.open).toBe(true)
    heatmap.click()
    await tick()
    expect(filters.show_heatmap).toBe(true)

    const column_checkboxes =
      column_menu.querySelectorAll<HTMLInputElement>(`.toggle-label input`)
    const checkbox_labels = [...column_menu.querySelectorAll(`.toggle-label`)].map(
      (label) => label.textContent?.trim(),
    )
    expect(checkbox_labels).toStrictEqual(sample_columns.map((column) => column.label))

    const [first_checkbox] = column_checkboxes
    expect(first_checkbox.checked).toBe(true)
    first_checkbox.click()
    expect(first_checkbox.checked).toBe(false)
    first_checkbox.click()
    expect(first_checkbox.checked).toBe(true)

    toggle_btn.click()
    expect(details.open).toBe(false)
  })
})
