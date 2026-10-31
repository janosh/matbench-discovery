import { parse_dependency_spec } from '#lib/environment.js'
import fs from 'node:fs'
import os from 'node:os'
import zlib from 'node:zlib'
import { build } from 'vite'
import { heading_anchor_html } from 'svelte-widgets/heading-anchors'
import { describe, expect, it, vi } from 'vite-plus/test'
import pkg from '../../package.json' with { type: 'json' }
import vite_config, { svelte_config } from '../../vite.config.js'

it(`bundles compressed JSON imports in production workers`, async () => {
  const root = fs.mkdtempSync(`${os.tmpdir()}/mbd-worker-build-`)
  try {
    fs.writeFileSync(`${root}/data.json.gz`, zlib.gzipSync(`{"symbol":"Cu"}`))
    fs.writeFileSync(
      `${root}/worker.js`,
      `import data from './data.json.gz'; postMessage(data.symbol)`,
    )
    fs.writeFileSync(
      `${root}/main.js`,
      `new Worker(new URL('./worker.js', import.meta.url), { type: 'module' })`,
    )
    await build({
      configFile: false,
      root,
      logLevel: `silent`,
      worker: vite_config.worker,
      build: { rolldownOptions: { input: `${root}/main.js` } },
    })
    const assets_dir = `${root}/dist/assets`
    const worker_file = fs
      .readdirSync(assets_dir)
      .find((file) => file.startsWith(`worker-`))
    if (!worker_file) throw new Error(`No worker chunk in ${assets_dir}`)
    const worker_code = fs.readFileSync(`${assets_dir}/${worker_file}`, `utf8`)
    // the gzipped payload is inlined as parsed JSON, not left as a binary import
    expect(worker_code).toContain(`JSON.parse(\`{"symbol":"Cu"}\`)`)
  } finally {
    fs.rmSync(root, { recursive: true, force: true })
  }
})

describe(`parse_dependency_spec`, () => {
  it.each([
    {
      // trim + extras + local version label; PyPI hrefs drop extras and the label
      dep: `  fairchem-core[torch-extras]==1.10.0+cu128  `,
      expected: {
        name: `fairchem-core[torch-extras]`,
        name_href: `https://pypi.org/project/fairchem-core/`,
        detail: `==1.10.0+cu128`,
        href: `https://pypi.org/project/fairchem-core/1.10.0`,
      },
    },
    {
      // non-== constraint (with spaces) links to project root
      dep: `mace-torch >= 0.3.16`,
      expected: {
        name: `mace-torch`,
        name_href: `https://pypi.org/project/mace-torch/`,
        detail: `>=0.3.16`,
        href: `https://pypi.org/project/mace-torch/`,
      },
    },
    {
      // git+ locator strips the git+ prefix; keeps URL fragment
      dep: `hienet @ git+https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
      expected: {
        name: `hienet`,
        name_href: `https://pypi.org/project/hienet/`,
        detail: `git+https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
        href: `https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet`,
      },
    },
    {
      // GitHub commit locators use the repository's browser URL shape
      dep: `tace @ git+https://github.com/xvzemin/tace@81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
      expected: {
        name: `tace`,
        name_href: `https://github.com/xvzemin/tace`,
        detail: `git+https://github.com/xvzemin/tace@81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
        href: `https://github.com/xvzemin/tace/commit/81f65a4c188bd09cec8d1419388f7afdcc1b6fd0`,
      },
    },
    {
      // symbolic GitHub revisions link to the corresponding repository tree
      dep: `aviary @ git+https://github.com/CompRhys/aviary.git@v0.1.0`,
      expected: {
        name: `aviary`,
        name_href: `https://github.com/CompRhys/aviary`,
        detail: `git+https://github.com/CompRhys/aviary.git@v0.1.0`,
        href: `https://github.com/CompRhys/aviary/tree/v0.1.0`,
      },
    },
    {
      // plain https locator used as-is
      dep: `aviary @ https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
      expected: {
        name: `aviary`,
        name_href: `https://pypi.org/project/aviary/`,
        detail: `https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
        href: `https://github.com/CompRhys/aviary/releases/tag/v0.1.0`,
      },
    },
    {
      // non-http locator falls back to PyPI
      dep: `pkg @ file:///tmp/wheel.whl`,
      expected: {
        name: `pkg`,
        name_href: `https://pypi.org/project/pkg/`,
        detail: `file:///tmp/wheel.whl`,
        href: `https://pypi.org/project/pkg/`,
      },
    },
    {
      dep: `ase`,
      expected: {
        name: `ase`,
        name_href: `https://pypi.org/project/ase/`,
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
  const { svelte_config: fresh_config } = await import(`../../vite.config.js`)
  expect(fresh_config.version).toStrictEqual(expected)
})

it(`first svelte preprocessor rewrites pkg.homepage links to site-internal paths`, async () => {
  const [strip_homepage] = svelte_config.preprocess
  const content = `<a href="${pkg.homepage}/models">models</a> ${pkg.homepage}`
  const result = await strip_homepage.markup?.({ content, filename: `readme.md` })
  expect(result?.code).toBe(`<a href="/models">models</a> `)
})

it(`manuscript preprocessing numbers figure labels and resolves forward references`, async () => {
  // find by name: Kit 3's sveltekit() appends its own warning preprocessor to the array
  const figure_markup = svelte_config.preprocess.find(
    (preprocessor) => preprocessor.name === `figure-refs`,
  )?.markup
  const filename = `site/src/routes/benchmarks/discovery/tmi/discovery-metric-figs.md`
  const content = fs.readFileSync(`../${filename}`, `utf8`)
  const result = await figure_markup?.({ content, filename })
  const labels = [...content.matchAll(/@label:(?<id>fig:[^\s]+)/g)]
  expect(labels).toHaveLength(5)
  for (const [idx, [, label]] of labels.entries()) {
    // self-contained anchor icon whose fragment keeps ':' literal so SvelteKit's
    // prerender crawler matches the id
    const anchor = heading_anchor_html(label)
    expect(anchor).toContain(`href="#${label}"`)
    expect(result?.code).toContain(
      `<strong id='${label}'>${anchor}Fig. ${idx + 1}</strong>`,
    )
  }
  expect(result?.code).not.toContain(`<use`)
  expect(result?.code).toContain(`<a href="#fig:cumulative-precision-recall">fig. 2</a>`)
  expect(result?.code).not.toMatch(/@(?:label|fig):/)
  expect(await figure_markup?.({ content, filename: `unrelated.md` })).toEqual({
    code: content,
  })
})
