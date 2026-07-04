import { goto } from '$app/navigation'
import type { GitHubActivityData } from '$lib/types'
import GitHubActivityScatter from '$lib/plot/GitHubActivityScatter.svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { get_scatter_plot_props, mount } from '../index'

const plot_mocks = vi.hoisted(() => ({
  ScatterPlot: vi.fn(),
}))

vi.mock(`matterviz`, () => ({ format_num: vi.fn(String) }))
vi.mock(`matterviz/plot`, () => ({ ScatterPlot: plot_mocks.ScatterPlot }))

type ScatterPlotProps = {
  series: {
    x: number[]
    y: number[]
    metadata: { name: string; repo: string; model_key?: string }[]
    color_values: number[]
    size_values: number[]
    point_label: { text: string; font_size: string; auto_placement: boolean }[]
  }[]
  x_axis: { label: string; format: string; range: [number, null] }
  y_axis: { label: string; format: string; range: [number, null] }
  color_bar: { title: string; tick_format: string }
  color_scale: { type: string }
  point_events: {
    onclick: (event: {
      point: { metadata?: { model_key?: string; name?: string } }
    }) => void
  }
  style?: string
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
  ...overrides,
})

describe(`GitHubActivityScatter`, () => {
  beforeEach(() => vi.clearAllMocks())

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
    expect(props.series).toHaveLength(1)
    expect(props.series[0]).toStrictEqual({
      x: [0, 100_000],
      y: [10, 999_999],
      markers: `points`,
      point_style: { fill: `#4dabf7` },
      metadata: [
        { name: `Zero Commits`, repo: `org/zero`, model_key: `zero-model` },
        { name: `Popular`, repo: `org/popular`, model_key: undefined },
      ],
      color_values: [1, 10_000],
      size_values: [1, 1000],
      point_label: [
        { text: `Zero Commits`, font_size: `11px`, auto_placement: true },
        { text: `Popular`, font_size: `11px`, auto_placement: true },
      ],
    })
    expect(props).toMatchObject({
      x_axis: { label: `GitHub Forks`, format: `,.0f`, range: [0, null] },
      y_axis: { label: `GitHub Stars`, format: `,.0f`, range: [0, null] },
      color_bar: { title: `Commits Last Year`, tick_format: `~s` },
      color_scale: { type: `log` },
    })
  })

  it(`passes rest props and disables labels when requested`, () => {
    mount(GitHubActivityScatter, {
      target: document.body,
      props: {
        github_data: [create_mock_github_data()],
        show_model_labels: false,
        style: `height: 800px;`,
      },
    })

    expect(get_scatter_plot_props(plot_mocks.ScatterPlot)).toMatchObject({
      series: [{ point_label: [] }],
      style: `height: 800px;`,
    })
  })

  it.each([
    {
      name: `model key`,
      overrides: { model_key: `preferred-slug`, name: `Display Name` },
      expected_path: `/models/preferred-slug`,
    },
    {
      name: `name fallback`,
      overrides: { name: `Name With Space` },
      expected_path: `/models/Name%20With%20Space`,
    },
  ])(`navigates to model page using $name`, ({ overrides, expected_path }) => {
    mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data(overrides)] },
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
    expect(goto).toHaveBeenCalledWith(expected_path)
  })
})
