// Per-model harmonic phonon "modes" assets for the kappa task page: the seekpath band
// path per phononDB material, converted into matterviz PhononModeData for the
// PhononModeExplorer. Only models whose kappa run stored the harmonic sidecar
// (`phonon_file` in the model YAML) get one. site/scripts/generate-kappa-parity-assets.py
// validates the shapes, so the client trusts the JSON without rechecking.
import type { Matrix3x3, Vec3 } from 'matterviz/math'
import type { Complex, PhononModeData, PhononPathSegment } from 'matterviz/spectral'
import { load_json_asset, parity_asset_resolver } from '../asset-loader'
import { kappa_parity_manifest } from './kappa-parity'

export interface KappaModesMaterial {
  primitive: {
    lattice: Matrix3x3
    symbols: string[]
    masses: number[]
    frac_coords: Vec3[]
  }
  segments: PhononPathSegment[]
  q_points: Vec3[]
  distances: number[]
  frequencies: number[][] // (n_q, n_bands) in THz
  // keyed by q-point index (JSON object keys are strings): [band][atom][xyz] -> [re, im]
  eigenvectors: Record<string, Complex[][][] | undefined>
}

export interface KappaModesModel {
  model_key: string
  materials: Record<string, KappaModesMaterial | undefined>
}

const {
  asset_url: kappa_modes_asset_url,
  model_asset: kappa_modes_asset,
  has_model: has_kappa_modes,
} = parity_asset_resolver(
  `kappa modes`,
  kappa_parity_manifest,
  import.meta.env.VITE_KAPPA_PARITY_ASSET_BASE_URL,
  `modes`,
)
export { has_kappa_modes }

export async function load_kappa_modes(model_key: string): Promise<KappaModesModel> {
  const { model } = await load_json_asset<{ model?: KappaModesModel }>(
    kappa_modes_asset_url(kappa_modes_asset(model_key)),
  )
  if (!model) throw new Error(`No kappa modes model ${model_key} in its asset`)
  if (model.model_key !== model_key) {
    throw new Error(
      `Invalid kappa modes model: expected ${model_key}, got ${model.model_key}`,
    )
  }
  return model
}

// q-points without stored eigenvectors get `eigenvector: null`; the explorer snaps to
// the nearest q-point that has one
export const build_phonon_mode_data = ({
  primitive,
  segments,
  q_points,
  distances,
  frequencies,
  eigenvectors,
}: KappaModesMaterial): PhononModeData => ({
  n_atoms: primitive.symbols.length,
  atoms: primitive.symbols.map((symbol, atom_idx) => ({
    symbol,
    mass: primitive.masses[atom_idx],
    coordinates: primitive.frac_coords[atom_idx],
  })),
  lattice: primitive.lattice,
  qpoints: q_points.map((q_position, q_idx) => ({
    q_position,
    distance: distances[q_idx],
    modes: frequencies[q_idx].map((frequency, band_idx) => ({
      frequency,
      eigenvector: eigenvectors[q_idx]?.[band_idx] ?? null,
    })),
  })),
  path_segments: segments,
})
