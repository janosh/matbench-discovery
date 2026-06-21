import { config } from '@janosh/vite-config'
import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import { load as yaml_load } from 'js-yaml'
import type { JSONSchema4 } from 'json-schema'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { execFileSync } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'
import zlib from 'node:zlib'
import type { Plugin } from 'vite'
import { defineConfig } from 'vite-plus'

// Load committed figure payloads (site/src/figs) as parsed ES modules, typed per payload
// in src/figs/payloads.d.ts. Two formats: <name>.json.gz (static, gzipped) and
// <name>.jsonl (multi-model, one model per line + a {"_base": {...}} shared-fields line,
// reassembled to { ...shared, models: [...] } - the format that git-merges cleanly).
// Both embed as JSON.parse('...') (cheap V8 parse; @__PURE__ drops unused payloads).
// .jsonl goes through attach_style ($lib/fig-helpers) so pages import pre-styled models.
const figure_payload_plugin = (): Plugin => ({
  name: `figure-payload`,
  load(id) {
    const file = id.split(`?`)[0]
    if (file.endsWith(`.json.gz`)) {
      const json = zlib.gunzipSync(fs.readFileSync(file)).toString(`utf8`)
      const code = `export default /* @__PURE__ */ JSON.parse(${JSON.stringify(json)})`
      return { code, map: null }
    }
    if (file.endsWith(`.jsonl`)) {
      const base: Record<string, unknown> = {}
      const models: unknown[] = []
      // trim each line (CRLF checkouts leave a trailing \r) before JSON.parse, matching
      // read_jsonl_payload's line.strip() in matbench_discovery/figs.py
      for (const raw of fs.readFileSync(file, `utf8`).split(`\n`)) {
        const line = raw.trim()
        if (!line) continue
        const entry = JSON.parse(line) as Record<string, unknown>
        if (`_base` in entry && Object.keys(entry).length === 1) {
          // oxlint-disable-next-line no-underscore-dangle -- _base is the shared-fields sentinel
          Object.assign(base, entry._base as Record<string, unknown>)
        } else models.push(entry)
      }
      const json = JSON.stringify({ ...base, models })
      const code = `import { attach_style } from '$lib/fig-helpers';
export default /* @__PURE__ */ attach_style(JSON.parse(${JSON.stringify(json)}))`
      return { code, map: null }
    }
    return null
  },
})

// Custom Vite plugin that watches for changes to *-schema.yml files and
// Automatically converts them to *-schema.d.ts files
function yaml_schema_to_typescript_plugin(): Plugin {
  // Convert *-schema.yml to *-schema.d.ts
  async function yaml_schema_to_ts(file: string): Promise<boolean> {
    try {
      const yaml_content = fs.readFileSync(file, `utf-8`)
      const file_dir = path.dirname(file)

      // Replace relative file paths in $refs with absolute file URIs
      const parsed_yaml = yaml_load(yaml_content) as JSONSchema4
      const base_name = path.basename(file, `.yml`)

      // Convert model schema to TypeScript
      const ts_name = {
        'model-schema': `ModelMetadata`,
        'dataset-schema': `DatasetRecord`,
        'label-schema': `Label`,
      }[base_name]

      const model_metadata_ts = await json_to_ts(
        parsed_yaml,
        ts_name ?? `missing ${base_name} schema`,
        {
          // Should match fmt config in vite.config.ts
          style: { semi: false, singleQuote: true, printWidth: 90 },
          // no-redundant-type-constituents fires on `string | 'missing'` unions from
          // literal-or-string schema fields and has no oxlint autofix, so suppress it
          // file-wide (index signatures are instead rewritten to Record<> via lint --fix below)
          bannerComment: `// This file is auto-generated from ${base_name}.yml. Do not edit directly.
// oxlint-disable typescript/no-redundant-type-constituents`,
          cwd: file_dir, // Important for $ref resolution between *-schema.yml files
          unreachableDefinitions: true,
        },
      )

      // Write the TypeScript interface file for model schema
      const dts_file = path.resolve(`./src/lib/${base_name}.d.ts`)
      fs.writeFileSync(dts_file, model_metadata_ts)

      // Rewrite json-schema-to-typescript index signatures (`{ [k: string]: T }`)
      // into `Record<string, T>` via the consistent-indexed-object-style autofix,
      // then format the generated schema file
      const vp_cmd = path.resolve(`./node_modules/.bin/vp`)
      execFileSync(vp_cmd, [`lint`, `--fix`, dts_file])
      execFileSync(vp_cmd, [`fmt`, `--write`, dts_file])
      return true
    } catch (error) {
      console.error(`Error processing schema file ${file}:`, error)
      return false
    }
  }

  return {
    name: `yaml-schema-to-typescript`,
    configureServer(server) {
      const schema_files = new Set(
        fs
          .readdirSync(`../tests`)
          .filter((file_name) => file_name.endsWith(`-schema.yml`))
          .map((file_name) => {
            const schema_file = path.resolve(`../tests/${file_name}`)
            void yaml_schema_to_ts(schema_file)
            server.watcher.add(schema_file) // Watch each schema file for changes
            return schema_file
          }),
      )

      server.watcher.on(`change`, (changed_file) => {
        if (schema_files.has(changed_file)) void yaml_schema_to_ts(changed_file)
      })
    },
  }
}

export default defineConfig({
  ...config, // shared lint/fmt/build from @janosh/vite-config (dotfiles)
  fmt: {
    ...config.fmt,
    ignorePatterns: [`src/routes/**/*.json`],
  },
  plugins: [
    sveltekit(),
    yaml_plugin({ extensions: [`.yml`, `.yaml`, `.cff`] }),
    yaml_schema_to_typescript_plugin(),
    figure_payload_plugin(),
  ],

  server: {
    fs: { allow: [`../..`] }, // Needed to import from $root
    port: 3000,
  },

  preview: {
    port: 3000,
  },

  test: {
    environment: `happy-dom`, // Faster than jsdom
    css: true,
    pool: `threads`,
    coverage: {
      reporter: [`text`, `json-summary`],
    },
    setupFiles: `tests/index.ts`,
    dir: `tests`,
    include: [`**/*.test.ts`, `**/*.test.svelte.ts`],
    // matterviz imports @spglib/moyo-wasm whose wasm binary can't load under vitest/node
    // (the .wasm references a missing _bg.js glue), so stub both the bare package and its
    // `/moyo_wasm_bg.wasm?url` subpath with a no-op mock
    alias: [
      {
        find: /^@spglib\/moyo-wasm(?:$|\/moyo_wasm_bg\.wasm)/,
        replacement: path.resolve(`tests/mocks/moyo-wasm.ts`),
      },
    ],
  },

  resolve: {
    conditions: process.env.VITEST ? [`browser`] : undefined,
  },
})
