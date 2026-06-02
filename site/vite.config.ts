import { config } from '@janosh/vite-config'
import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import { load as yaml_load } from 'js-yaml'
import type { JSONSchema4 } from 'json-schema'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { execFileSync } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'
import type { Plugin } from 'vite'
import { defineConfig } from 'vite-plus'

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
    ignorePatterns: [`src/figs/**/*.svelte`, `src/routes/**/*.json`],
  },
  // Shared rules/plugins/categories live in @janosh/vite-config (dotfiles).
  // Append only matbench-discovery-specific ignore dirs and rule overrides here.
  lint: {
    ...config.lint,
    ignorePatterns: [...config.lint.ignorePatterns, `src/figs/**`],
  },
  staged: {
    '*': `codespell --ignore-words-list falsy --check-filenames`,
    '*.test.ts': `sh -c '! grep -E "(test|describe)\\.only\\(" "$@"' --`,
    '*.{js,ts,svelte,html,css,md,json,yaml}': `vp check --fix`,
    '*.{ts,svelte}': `sh -c 'pnpm exec svelte-kit sync && pnpm exec svelte-check-rs --threshold error'`,
  },
  plugins: [
    sveltekit(),
    yaml_plugin({ extensions: [`.yml`, `.yaml`, `.cff`] }),
    yaml_schema_to_typescript_plugin(),
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
    pool: `threads`, // Parallel test execution
    coverage: {
      reporter: [`text`, `json-summary`],
    },
    setupFiles: `tests/index.ts`,
    dir: `tests`,
    include: [`**/*.test.ts`, `**/*.test.svelte.ts`],
  },

  resolve: {
    conditions: process.env.VITEST ? [`browser`] : undefined,
    alias: process.env.VITEST
      ? [
          // Mock wasm-dependent modules to avoid loading issues in jsdom
          {
            find: /^@spglib\/moyo-wasm.*/,
            replacement: new URL(`tests/mocks/moyo-wasm.ts`, import.meta.url).pathname,
          },
        ]
      : [],
  },
})
