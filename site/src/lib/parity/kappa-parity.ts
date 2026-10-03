import type kappa_analysis_data from '#figs/kappa-103-analysis.jsonl'
import type { PhononDos } from 'matterviz/spectral'
import {
  assert_array_length,
  load_json_asset,
  load_parity_model,
  parity_asset_resolver,
} from '../asset-loader.js'
import type { ParityBase, ParityModel, ParityPoint } from '../asset-loader.js'
import kappa_parity_manifest_json from './kappa-parity-manifest.json'
import { is_finite_num } from '../metrics.js'

export type KappaAnalysis = typeof kappa_analysis_data

let analysis_promise: Promise<KappaAnalysis> | undefined

// Lazily load and cache the shared κ-103 per-material analysis payload.
export const load_kappa_analysis = (): Promise<KappaAnalysis> =>
  (analysis_promise ??= import(`#figs/kappa-103-analysis.jsonl`).then(
    (module) => module.default,
  ))

// Map material IDs to one model's per-material κ_SRME values.
export async function load_kappa_srme_map(
  model_key: string,
): Promise<Map<string, number | null> | undefined> {
  const analysis = await load_kappa_analysis()
  const model_analysis = analysis.models.find(
    (model_data) => model_data.model_key === model_key,
  )
  if (!model_analysis) return undefined
  return new Map(
    analysis.material_ids.map((material_id, idx) => [
      material_id,
      model_analysis.srme[idx],
    ]),
  )
}

export const kappa_parity_manifest = kappa_parity_manifest_json
// raw phonon DOS as stored in assets (histogram of mesh frequencies in THz)
interface RawDos {
  frequencies: number[]
  densities: number[]
}

export interface KappaParityBase extends ParityBase {
  kappa_dft: (number | null)[]
  n_sites: (number | null)[]
  spacegroups: (number | null)[]
  dft_dos: Record<string, RawDos | undefined>
}

export interface KappaParityModel extends ParityModel {
  kappa_ml: (number | null)[]
  ml_dos: Record<string, RawDos | undefined>
}

// type alias (not interface) so it satisfies ScatterPlot's Record<string, unknown>
// metadata constraint without casts
export type KappaParityPoint = ParityPoint & {
  kappa_dft: number
  kappa_ml: number
  sre: number // symmetric relative error of scalar conductivity
  n_sites: number | null
  spacegroup: number | null
}

interface KappaParitySeries {
  x: number[]
  y: number[]
  points: KappaParityPoint[]
}

export const {
  asset_url: kappa_parity_asset_url,
  model_asset: kappa_model_asset,
  has_model: has_kappa_parity_model,
} = parity_asset_resolver(
  `kappa`,
  kappa_parity_manifest,
  import.meta.env.VITE_KAPPA_PARITY_ASSET_BASE_URL,
)

export const load_kappa_parity_base = (): Promise<KappaParityBase> =>
  load_json_asset<KappaParityBase>(
    kappa_parity_asset_url(kappa_parity_manifest.base.asset),
    (base) => {
      for (const key of [
        `material_ids`,
        `formulas`,
        `kappa_dft`,
        `n_sites`,
        `spacegroups`,
      ] as const) {
        assert_array_length(
          `kappa parity ${key}`,
          base[key],
          kappa_parity_manifest.row_count,
        )
      }
    },
  )

export const load_kappa_parity_model = (model_key: string): Promise<KappaParityModel> =>
  load_parity_model<KappaParityModel>(
    `kappa`,
    kappa_parity_asset_url(kappa_model_asset(model_key)),
    model_key,
    `kappa_ml`,
    kappa_parity_manifest.row_count,
  )

export function get_kappa_parity_point(
  base: KappaParityBase,
  model: KappaParityModel,
  row_idx: number,
): KappaParityPoint | null {
  const x = base.kappa_dft[row_idx]
  const y = model.kappa_ml[row_idx]
  // require positive conductivities: physically meaningful and log-scale safe
  if (!is_finite_num(x) || !is_finite_num(y) || x <= 0 || y <= 0) return null

  return {
    material_id: base.material_ids[row_idx] ?? `unknown-${row_idx}`,
    formula: base.formulas[row_idx] ?? ``,
    kappa_dft: x,
    kappa_ml: y,
    sre: Math.abs((2 * (y - x)) / (x + y)),
    n_sites: base.n_sites[row_idx] ?? null,
    spacegroup: base.spacegroups[row_idx] ?? null,
  }
}

export function build_kappa_parity_series(
  base: KappaParityBase,
  model: KappaParityModel,
): KappaParitySeries {
  const points = base.material_ids
    .map((_material_id, row_idx) => get_kappa_parity_point(base, model, row_idx))
    .filter((point): point is KappaParityPoint => point !== null)
  const x = points.map((point) => point.kappa_dft)
  const y = points.map((point) => point.kappa_ml)
  return { x, y, points }
}

export function as_phonon_dos(dos: RawDos | undefined): PhononDos | null {
  if (!dos?.frequencies?.length || !dos.densities?.length) return null
  if (dos.frequencies.length !== dos.densities.length) return null
  return { type: `phonon`, frequencies: dos.frequencies, densities: dos.densities }
}

// Whether the DOS carries weight at negative frequencies, i.e. the structure is
// dynamically unstable under this model. The thermal plot deliberately keeps those modes
// in its DOS integral so the curves come out depressed; this lets the UI say why,
// instead of showing an unexplained sub-3k_B heat capacity
export const has_imaginary_modes = (dos: PhononDos): boolean =>
  dos.frequencies.some((freq, idx) => freq < 0 && dos.densities[idx] > 0)
