// Type declarations for Node.js modules
/// <reference types="node" />

import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import yaml from 'js-yaml'
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
      const parsed_yaml = yaml.load(yaml_content) as JSONSchema4
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
          bannerComment: `// This file is auto-generated from ${base_name}.yml. Do not edit directly.`,
          cwd: file_dir, // Important for $ref resolution between *-schema.yml files
          unreachableDefinitions: true,
        },
      )

      // Write the TypeScript interface file for model schema
      const dts_file = path.resolve(`./src/lib/${base_name}.d.ts`)
      fs.writeFileSync(dts_file, model_metadata_ts)

      // Format generated schema file
      const vp_cmd = path.resolve(`./node_modules/.bin/vp`)
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
  fmt: {
    printWidth: 90,
    semi: false,
    singleQuote: true,
    ignorePatterns: [
      `src/figs/**/*.svelte`,
      `src/figs/**/*.json`,
      `src/routes/**/*.json`,
      `tests/**`, // oxfmt lowercases first char of describe/it/test strings
    ],
  },
  lint: {
    categories: {
      correctness: `error`,
      suspicious: `error`,
      pedantic: `error`,
      perf: `error`,
    },
    ignorePatterns: [
      `build/`,
      `.svelte-kit/`,
      `src/figs/`,
      `dist/`,
      `src/lib/*.d.ts`,
      `scripts/`,
    ],
    options: {
      typeAware: true,
      typeCheck: true,
    },
    plugins: [`oxc`, `typescript`, `unicorn`, `import`, `jest`],
    rules: {
      'no-console': [`error`, { allow: [`warn`, `error`] }],
      'no-unused-vars': `off`,
      '@typescript-eslint/no-unused-vars': [
        `error`,
        { argsIgnorePattern: `^_`, varsIgnorePattern: `^_` },
      ],
      'no-self-assign': `off`,
      'no-await-in-loop': `off`,
      'no-shadow': `off`,
      'prefer-const': `off`,
      'no-unassigned-vars': `off`,
      '@typescript-eslint/no-unnecessary-condition': `off`,
      '@typescript-eslint/consistent-type-imports': `off`,
      'eslint-plugin-unicorn/consistent-function-scoping': `off`,
      'eslint-plugin-unicorn/no-new-array': `off`,
      'eslint-plugin-unicorn/prefer-array-find': `off`,
      'eslint-plugin-import/no-self-import': `off`,
      'eslint-plugin-import/no-unassigned-import': `off`,
      'eslint-plugin-import/no-named-as-default-member': `off`,
      '@typescript-eslint/no-unsafe-argument': `off`,
      '@typescript-eslint/no-unsafe-assignment': `off`,
      '@typescript-eslint/no-unsafe-call': `off`,
      '@typescript-eslint/no-unsafe-member-access': `off`,
      '@typescript-eslint/no-unsafe-return': `off`,
      'no-confusing-void-expression': `off`,
      'max-lines-per-function': `off`,
      'strict-boolean-expressions': `off`,
      'eslint-plugin-jest/no-conditional-in-test': `off`,
      'eslint-plugin-jest/no-conditional-expect': `off`,
      'eslint-plugin-jest/no-disabled-tests': `off`,
      'eslint-plugin-jest/valid-expect': `off`,
      '@typescript-eslint/no-unsafe-type-assertion': `off`,
      '@typescript-eslint/only-throw-error': `off`,
      '@typescript-eslint/no-misused-promises': `off`,
      '@typescript-eslint/no-deprecated': `off`,
      '@typescript-eslint/await-thenable': `off`,
      'eslint-plugin-unicorn/no-array-callback-reference': `off`,
      'eslint-plugin-unicorn/no-useless-undefined': `off`,
      'eslint-plugin-unicorn/no-object-as-default-parameter': `off`,
      'eslint-plugin-import/max-dependencies': `off`,
      'no-inline-comments': `off`,
      'no-negated-condition': `off`,
      'no-warning-comments': `off`,
      'no-promise-executor-return': `off`,
      'oxc/no-map-spread': `off`,
      'max-lines': `off`,
      'max-depth': `off`,
    },
  },
  staged: {
    '*': `codespell --ignore-words-list falsy --check-filenames`,
    '*.test.ts': `sh -c '! grep -E "(test|describe)\\.only\\(" "$@"' --`,
    '*.{js,ts,svelte,html,css,md,json,yaml}': `vp check --fix`,
    '*.{ts,svelte}': `sh -c 'npx svelte-kit sync && npx svelte-check-rs --threshold error'`,
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
    server: {
      deps: {
        // Force Vitest to inline packages with directory imports to handle them properly
        inline: [`matterviz`, `svelte-multiselect`, `@threlte/core`, `@threlte/extras`],
      },
    },
  },

  resolve: {
    conditions: process.env.VITEST ? [`browser`] : undefined,
    alias: process.env.VITEST
      ? [
          // Mock wasm-dependent modules to avoid loading issues in jsdom
          {
            find: /^@spglib\/moyo-wasm.*/,
            replacement: new URL(`./tests/mocks/moyo-wasm.ts`, import.meta.url).pathname,
          },
        ]
      : [],
  },
})
