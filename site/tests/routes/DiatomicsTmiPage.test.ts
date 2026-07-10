import DiatomicsTmiPage from '$routes/tasks/diatomics/tmi/+page.svelte'
import { tick } from 'svelte'
import { afterEach, expect, it, vi } from 'vitest'
import { get_scatter_plot_props, mount_with_url } from '../index'

const plot_mocks = vi.hoisted(() => ({
  ScatterPlot: vi.fn(),
}))

vi.mock(`matterviz`, () => ({ ScatterPlot: plot_mocks.ScatterPlot }))

type PlotSeries = { id: string; x: number[]; y: number[] }

const plotted_series = (): PlotSeries[] =>
  (get_scatter_plot_props(plot_mocks.ScatterPlot) as { series: PlotSeries[] }).series

function button_for(text: string, selector = `button`): HTMLButtonElement {
  const button = [...document.querySelectorAll<HTMLButtonElement>(selector)].find(
    (candidate) => candidate.textContent?.trim() === text,
  )
  if (!button) throw new Error(`${text} button not found`)
  return button
}

afterEach(() => vi.unstubAllGlobals())

it(`filters TMI elements and functionals while preserving curve gaps`, async () => {
  vi.stubGlobal(`IntersectionObserver`, undefined)
  const magmom_curve = {
    distances: [1, 2],
    magmoms: [[1, -1] as [number, number], null],
    spin_candidates: [`afm`, null],
  }
  await mount_with_url(
    DiatomicsTmiPage,
    `http://localhost/tasks/diatomics/tmi?elements=halogen`,
    {
      props: {
        data: {
          magmom_curves: {
            'F-F': { PBE: magmom_curve, r2SCAN: magmom_curve },
          },
        },
      },
    },
  )

  expect(button_for(`Halogens`).getAttribute(`aria-pressed`)).toBe(`true`)
  const pbe = button_for(`PBE`, `.legend button`)
  expect(plotted_series()).toHaveLength(4)
  expect(plotted_series()[0]).toMatchObject({ x: [1, 2], y: [1, Number.NaN] })

  pbe.click()
  await tick()
  expect(plotted_series().map((series) => series.id)).toEqual([
    `F-F-r2SCAN-atom1`,
    `F-F-r2SCAN-atom2`,
  ])

  button_for(`Lanthanides`).click()
  await tick()
  expect(new URL(location.href).searchParams.get(`elements`)).toBe(`lanthanide`)
})
