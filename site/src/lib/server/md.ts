import type { ModelData } from '$lib/types'
import { readFile } from 'node:fs/promises'
import { isAbsolute, resolve } from 'node:path'
import { gunzipSync } from 'node:zlib'

const repo_root = resolve(import.meta.dirname, `../../../..`)

export type MdPerSystemRow = Record<string, number | string>

// Read a model's per-system MD metrics from its committed pred_file CSV (one row per
// DynaMat system). Returns null when the model has no MD data or the CSV is absent
// (e.g. a pred_file pending re-upload). Build-time only (node:fs).
export async function read_md_per_system(
  model: Pick<ModelData, `metrics`>,
  root_dir: string = repo_root,
): Promise<MdPerSystemRow[] | null> {
  const md = model.metrics?.md
  if (!md || typeof md !== `object` || !md.pred_file) return null

  const csv_path = isAbsolute(md.pred_file)
    ? md.pred_file
    : resolve(root_dir, md.pred_file)
  let csv: string
  try {
    const bytes = await readFile(csv_path)
    csv =
      bytes[0] === 0x1f && bytes[1] === 0x8b
        ? gunzipSync(bytes).toString(`utf-8`)
        : bytes.toString(`utf-8`)
  } catch (error) {
    if (error instanceof Error && `code` in error && error.code === `ENOENT`) return null
    throw error
  }

  const [header, ...lines] = csv.trim().split(`\n`)
  const cols = header.split(`,`)
  return lines.map((line) => {
    const row: MdPerSystemRow = {}
    for (const [idx, raw] of line.split(`,`).entries()) {
      if (raw === ``) continue // missing value (e.g. stress-less system pressure)
      const num = Number(raw)
      // per-system energy/force RMSEs are stored in eV; report meV to match the
      // model-level MD table units
      row[cols[idx]] = Number.isFinite(num)
        ? [`energy_rmse`, `force_rmse`].includes(cols[idx])
          ? num * 1000
          : num
        : raw
    }
    return row
  })
}
