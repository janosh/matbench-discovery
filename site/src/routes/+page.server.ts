import model_schema from '$root/tests/model-schema.yml'
import { prettier as prettier_config } from '$site/package.json'
import fs from 'fs'
import { compile as json_to_ts } from 'json-schema-to-typescript'

// keep model-schema.d.ts in sync with model-schema.yml (source of truth)
// i.e. use json-schema-to-typescript to auto-convert YAML schema to TypeScript interface
const model_metadata_ts = await json_to_ts(model_schema, `ModelMetadata`, {
  style: prettier_config,
})
// prettier format model_md_type
const dts_out_file = `src/lib/model-schema.d.ts`
fs.writeFileSync(dts_out_file, model_metadata_ts)
