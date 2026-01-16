import { GitHubActivityScatter } from '$lib'
import type { GitHubActivityData } from '$lib/types'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { doc_query } from '../index'

vi.mock(`matterviz`, () => ({
  format_num: (val: number) => String(val),
  Icon: vi.fn(),
}))
vi.mock(`matterviz/plot`, () => ({ ScatterPlot: vi.fn() }))

function create_mock_github_data(
  overrides: Partial<GitHubActivityData> = {},
): GitHubActivityData {
  return {
    name: `Test Model`,
    repo: `test-org/test-model`,
    stars: 1500,
    forks: 250,
    commits_last_year: 120,
    contributors: 15,
    ...overrides,
  }
}

describe(`GitHubActivityScatter`, () => {
  beforeEach(() => {
    document.exitFullscreen = vi.fn().mockResolvedValue(undefined)
    Element.prototype.requestFullscreen = vi.fn().mockResolvedValue(undefined)
  })

  it.each([
    [`empty data`, []],
    [`single point`, [create_mock_github_data()]],
    [`multiple points`, [
      create_mock_github_data({ name: `A`, stars: 5000 }),
      create_mock_github_data({ name: `B`, stars: 2500 }),
    ]],
  ])(`renders with %s`, (_, github_data) => {
    mount(GitHubActivityScatter, { target: document.body, props: { github_data } })
    expect(document.body.firstElementChild).toBeDefined()
  })

  it(`has bleed-1400 class when not fullscreen`, () => {
    mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data()] },
    })
    expect(doc_query<HTMLDivElement>(`div.bleed-1400`)).toBeDefined()
  })

  it.each([
    [`show_model_labels=false`, { show_model_labels: false }],
    [`custom style`, { style: `height: 800px;` }],
  ])(`handles %s prop`, (_, extra_props) => {
    const instance = mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data()], ...extra_props },
    })
    expect(instance).toBeDefined()
  })

  it(`filters out data points with null values`, () => {
    const instance = mount(GitHubActivityScatter, {
      target: document.body,
      props: {
        github_data: [
          create_mock_github_data(),
          {
            name: `Incomplete`,
            repo: `x`,
            stars: null,
            forks: 100,
            commits_last_year: 50,
            contributors: 5,
          } as unknown as GitHubActivityData,
        ],
      },
    })
    expect(instance).toBeDefined()
  })

  it.each([
    [`zero values`, { stars: 0, forks: 0, commits_last_year: 0, contributors: 1 }],
    [`large values`, {
      stars: 999999,
      forks: 100000,
      commits_last_year: 10000,
      contributors: 1000,
    }],
  ])(`handles edge case: %s`, (_, values) => {
    const instance = mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data(values)] },
    })
    expect(instance).toBeDefined()
  })

  it(`handles various repo name formats`, () => {
    const instance = mount(GitHubActivityScatter, {
      target: document.body,
      props: {
        github_data: [
          create_mock_github_data({ repo: `org/simple` }),
          create_mock_github_data({ repo: `org-name/with-dashes` }),
          create_mock_github_data({ repo: `OrgName/CamelCase` }),
        ],
      },
    })
    expect(instance).toBeDefined()
  })

  it(`responds to fullscreenchange event without error`, () => {
    mount(GitHubActivityScatter, {
      target: document.body,
      props: { github_data: [create_mock_github_data()] },
    })
    globalThis.dispatchEvent(new Event(`fullscreenchange`))
    expect(document.body.firstElementChild).toBeDefined()
  })
})
