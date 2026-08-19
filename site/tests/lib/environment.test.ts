import { parse_dependency_spec } from '$lib/environment'
import { describe, expect, it, vi } from 'vitest'

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

it(`uses SvelteKit's stable default version in production`, async () => {
  vi.stubEnv(`NODE_ENV`, `production`)
  vi.resetModules()
  try {
    const { default: svelte_config } = await import(`../../svelte.config`)
    expect(svelte_config.kit?.version).toBeUndefined()
  } finally {
    vi.unstubAllEnvs()
    vi.resetModules()
  }
})
