import { goto } from '$app/navigation'
import type { GitHubActivityData } from '#lib/types.js'
import GitHubActivityScatter from '#lib/plot/GitHubActivityScatter.svelte'
import { beforeEach, describe, expect, it, vi } from 'vite-plus/test'
import { get_scatter_plot_props, mount } from '../index.js'

// Svelte components and snippets both render as (anchor, ...prop getters), so the mock
// can render the tooltip snippet for one hovered point when a test sets tooltip_point.
const plot_mocks = vi.hoisted(() => {
  const state: { tooltip_point?: Record<string, unknown> } = {}
  type MockProps = { tooltip?: (anchor: Node, get_point: () => unknown) => void }
  return {
    state,
    ScatterPlot: vi.fn((anchor: Node, props: MockProps) => {
      const { tooltip_point } = state
      if (tooltip_point) props.tooltip?.(anchor, () => tooltip_point)
    }),
  }
})

vi.mock(`matterviz/labels`, () => ({ format_num: vi.fn(String) }))
vi.mock(`matterviz/plot`, () => ({ ScatterPlot: plot_mocks.ScatterPlot }))

type ScatterMetadata = { name?: string; model_key?: string }
type ScatterPlotProps = {
  series: { metadata: ScatterMetadata[] }[]
  x_axis: unknown
  point_events: { onclick: (event: { point: { metadata?: ScatterMetadata } }) => void }
}

const create_mock_github_data = (
  overrides: Partial<GitHubActivityData> = {},
): GitHubActivityData => ({
  name: `Test Model`,
  repo: `test-org/test-model`,
  stars: 1500,
  forks: 250,
  commits_last_year: 120,
  contributors: 15,
  model_key: `test-model`,
  ...overrides,
})

describe(`GitHubActivityScatter`, () => {
  beforeEach(() => {
    vi.clearAllMocks()
    plot_mocks.state.tooltip_point = undefined
  })

  it(`filters invalid activity rows and maps metrics to ScatterPlot props`, () => {
    mount(GitHubActivityScatter, {
      target: document.body,
      props: {
        github_data: [
          create_mock_github_data({
            name: `Zero Commits`,
            repo: `org/zero`,
            model_key: `zero-model`,
            stars: 10,
            forks: 0,
            commits_last_year: 0,
            contributors: 1,
          }),
          create_mock_github_data({
            name: `Popular`,
            repo: `org/popular`,
            stars: 999_999,
            forks: 100_000,
            commits_last_year: 10_000,
            contributors: 1000,
          }),
          create_mock_github_data({
            name: `Incomplete`,
            repo: `x`,
            stars: Number.NaN,
            forks: 100,
            commits_last_year: 50,
            contributors: 5,
          }),
          create_mock_github_data({
            name: `Infinite`,
            stars: Number.POSITIVE_INFINITY,
          }),
        ],
      },
    })

    expect(plot_mocks.ScatterPlot).toHaveBeenCalledTimes(1)
    const props = get_scatter_plot_props(plot_mocks.ScatterPlot) as ScatterPlotProps
    expect(props.series).toStrictEqual([
      {
        x: [0, 100_000],
        y: [10, 999_999],
        markers: `points`,
        point_style: { fill: `#4dabf7` },
        metadata: [
          {
            name: `Zero Commits`,
            repo: `org/zero`,
            model_key: `zero-model`,
            commits_last_year: 0,
            contributors: 1,
          },
          {
            name: `Popular`,
            repo: `org/popular`,
            model_key: `test-model`,
            commits_last_year: 10_000,
            contributors: 1000,
          },
        ],
        color_values: [1, 10_000],
        size_values: [1, 1000],
        point_label: [
          { text: `Zero Commits`, font_size: `11px`, auto_placement: true },
          { text: `Popular`, font_size: `11px`, auto_placement: true },
        ],
      },
    ])
    expect(props).toMatchObject({
      y_axis: { label: `GitHub Stars`, format: `,.0f`, range: [0, null] },
      color_bar: { title: `Commits Last Year`, tick_format: `~s` },
      color_scale: { type: `log` },
    })
    // exact x_axis: no label_shift, since an inward title shift puts the label over the data
    expect(props.x_axis).toStrictEqual({
      label: `GitHub Forks`,
      format: `,.0f`,
      range: [0, null],
    })
  })

  it(`passes rest props through and navigates to the clicked model's page`, () => {
    const overrides = { model_key: `slug with space`, name: `Display Name` }
    mount(GitHubActivityScatter, {
      target: document.body,
      props: {
        github_data: [create_mock_github_data(overrides)],
        style: `height: 800px;`,
      },
    })
    expect(get_scatter_plot_props(plot_mocks.ScatterPlot)).toMatchObject({
      style: `height: 800px;`,
    })

    const scatter_plot_props = get_scatter_plot_props(
      plot_mocks.ScatterPlot,
    ) as ScatterPlotProps
    const metadata = scatter_plot_props.series[0]?.metadata[0]
    expect(metadata).toMatchObject(overrides)
    scatter_plot_props.point_events.onclick({
      point: { metadata },
    })
    expect(goto).toHaveBeenCalledTimes(1)
    expect(goto).toHaveBeenCalledWith(`/models/slug%20with%20space`)
  })

  it(`shows commits and contributors from point metadata in the tooltip`, () => {
    plot_mocks.state.tooltip_point = {
      x_formatted: `250`,
      y_formatted: `1,500`,
      metadata: { name: `Test Model`, commits_last_year: 120, contributors: 15 },
    }
    mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data()] },
    })
    expect(document.body.textContent?.replaceAll(/\s+/g, ` `).trim()).toBe(
      `Test Model GitHub Forks: 250 GitHub Stars: 1,500 Commits Last Year: 120 Contributors: 15`,
    )
  })
})
