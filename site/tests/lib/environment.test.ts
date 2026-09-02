import { parse_dependency_spec } from '$lib/environment'
import { describe, expect, it, vi } from 'vitest'
import pkg from '../../package.json' with { type: 'json' }
import { svelte_config } from '../../vite.config'

describe(`parse_dependency_spec`, () => {
  it.each([
    {
      // trim + extras + local version label; PyPI href drops extras
      dep: `  fairchem-core[torch-extras]==1.10.0+cu128  `,
      expected: {
        name: `fairchem-core[torch-extras]`,
        detail: `==1.10.0+cu128`,
        href: `https://pypi.org/project/fairchem-core/1.10.0+cu128`,
      },
    },
    {
      // non-== constraint (with spaces) links to project root
      dep: `mace-torch >= 0.3.16`,
      expected: {
        name: `mace-torch`,
        detail: `>=0.3.16`,
        href: `https://pypi.org/project/mace-torch/`,
      },
    },
    {
      // git+ locator strips the git+ prefix; keeps URL fragment
      dep: `hienet @ git+https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
      expected: {
        name: `hienet`,
        detail: `git+https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
        href: `https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
      },
    },
    {
      // GitHub commit locators use the repository's browser URL shape
      dep: `tace @ git+https://github.com/xvzemin/tace@81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
      expected: {
        name: `tace`,
        detail: `git+https://github.com/xvzemin/tace@81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
        href: `https://github.com/xvzemin/tace/commit/81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
      },
    },
    {
      // symbolic GitHub revisions link to the corresponding repository tree
      dep: `aviary @ git+https://github.com/CompRhys/aviary.git@v0.1.0`,
      expected: {
        name: `aviary`,
        detail: `git+https://github.com/CompRhys/aviary.git@v0.1.0`,
        href: `https://github.com/CompRhys/aviary/tree/v0.1.0`,
      },
    },
    {
      // plain https locator used as-is
      dep: `aviary @ https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
      expected: {
        name: `aviary`,
        detail: `https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
        href: `https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
      },
    },
    {
      // non-http locator falls back to PyPI
      dep: `pkg @ file:///tmp/wheel.whl`,
      expected: {
        name: `pkg`,
        detail: `file:///tmp/wheel.whl`,
        href: `https://pypi.org/project/pkg/`,
      },
    },
    {
      dep: `ase`,
      expected: {
        name: `ase`,
        detail: ``,
        href: `https://pypi.org/project/ase/`,
      },
    },
  ])(`parses $dep`, ({ dep, expected }) => {
    expect(parse_dependency_spec(dep)).toStrictEqual(expected)
  })
})

// svelte_config is evaluated at module load, so re-import after stubbing NODE_ENV
it.each([
  { node_env: `production`, expected: undefined },
  { node_env: `development`, expected: { name: `dev` } },
])(`kit version is $expected when NODE_ENV=$node_env`, async ({ node_env, expected }) => {
  vi.stubEnv(`NODE_ENV`, node_env)
  vi.resetModules()
  try {
    const { svelte_config: fresh_config } = await import(`../../vite.config`)
    expect(fresh_config.version).toStrictEqual(expected)
  } finally {
    vi.unstubAllEnvs()
    vi.resetModules()
  }
})

it(`first svelte preprocessor rewrites pkg.homepage links to site-internal paths`, async () => {
  const [strip_homepage] = svelte_config.preprocess
  const content = `<a href="${pkg.homepage}/models">models</a> ${pkg.homepage}`
  const result = await strip_homepage.markup?.({ content, filename: `readme.md` })
  expect(result?.code).toBe(`<a href="/models">models</a> `)
})
