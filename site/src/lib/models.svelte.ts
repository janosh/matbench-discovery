import { DEFAULT_CPS_CONFIG, calculate_cps } from './metrics'
import type { ModelData } from './types'

// Make CPS_CONFIG reactive (using Svelte 5 runes)
export const CPS_CONFIG = $state({ ...DEFAULT_CPS_CONFIG })

export const MODEL_METADATA_PATHS = import.meta.glob(`$root/models/[^_]**/[^_]*.yml`, {
  eager: true,
  import: `default`,
}) as Record<string, ModelData>

// visually distinct color palette
export const MODEL_COLORS = [
  `#4285F4`, // blue
  `#EA4335`, // red
  `#FBBC05`, // yellow
  `#34A853`, // green
  `#8A2BE2`, // blueviolet
  `#FF7F50`, // coral
  `#1E90FF`, // dodgerblue
  `#FF1493`, // deeppink
  `#32CD32`, // limegreen
  `#FF8C00`, // darkorange
  `#9370DB`, // mediumpurple
  `#3CB371`, // mediumseagreen
  `#DC143C`, // crimson
  `#6495ED`, // cornflowerblue
  `#FFD700`, // gold
  `#8B008B`, // darkmagenta
  `#00CED1`, // darkturquoise
  `#FF4500`, // orangered
  `#2E8B57`, // seagreen
  `#BA55D3`, // mediumorchid
]

export const MODELS = $state(
  Object.entries(MODEL_METADATA_PATHS)
    .filter(
      // ignore models that aren't completed
      ([_key, metadata]) => (metadata?.status ?? `complete`) == `complete`,
    )
    .map(([key, metadata], index) => {
      // Assign color to each model for consistent coloring across plots
      const model_color = MODEL_COLORS[index % MODEL_COLORS.length]

      return {
        ...metadata,
        dirname: key.split(`/`)[2],
        metadata_file: key.replace(/^..\//, ``),
        color: model_color,
        CPS: NaN, // Initial CPS placeholder
      }
    }),
) as ModelData[]

// Update CPSs of models based on current CPS weights
export function update_models_cps() {
  MODELS.forEach((model: ModelData) => {
    // Extract required metrics for CPS calculation
    const f1 = model.metrics?.discovery?.[`unique_prototypes`]?.F1
    const rmsd =
      model.metrics?.geo_opt && typeof model.metrics.geo_opt !== `string`
        ? model.metrics.geo_opt[`symprec=1e-5`]?.rmsd
        : undefined
    const kappa =
      model.metrics?.phonons && typeof model.metrics.phonons !== `string`
        ? model.metrics.phonons.kappa_103?.κ_SRME !== undefined
          ? Number(model.metrics.phonons.kappa_103.κ_SRME)
          : undefined
        : undefined

    // Calculate and update CPS
    model.CPS = calculate_cps(f1, rmsd, kappa, CPS_CONFIG) ?? NaN
  })
}

// Calculate initial CPS for all models
update_models_cps()
