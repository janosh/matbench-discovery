import dft_references from '$lib/diatomics-dft.json.gz'
import type { PageServerLoad } from './$types'

// functional -> formula -> per-point magmom data (subset of the bundled DFT reference)
type DftMagmomCurve = {
  distances: number[]
  magmoms: ([number, number] | null)[]
  spin_candidates: (string | null)[]
}

const DFT_REFERENCES = dft_references as Record<string, Record<string, DftMagmomCurve>>

export const load: PageServerLoad = () => {
  // per-atom site-projected moments vs separation for every element/functional. This
  // is the spin-state debugging view Andrew Rosen suggested: discontinuities in
  // atom-wise moments flag SCF spin-state hops the total energy curve can hide.
  const magmom_curves: Record<string, Record<string, DftMagmomCurve>> = {}
  for (const [functional, formulas] of Object.entries(DFT_REFERENCES)) {
    for (const [formula, curve] of Object.entries(formulas)) {
      if (!curve.magmoms.some((moments) => moments !== null)) continue
      magmom_curves[formula] ??= {}
      magmom_curves[formula][functional] = {
        distances: curve.distances,
        magmoms: curve.magmoms,
        spin_candidates: curve.spin_candidates,
      }
    }
  }
  return { magmom_curves }
}
