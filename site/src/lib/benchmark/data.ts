import DATASETS from '$data/datasets.yml'
import data_files from '$pkg/data-files.yml'
import pkg from '$site/package.json'
import { format_num } from 'matterviz/labels'
import {
  Search,
  Download,
  Info,
  Molecule,
  Phonons,
  RulerSquareCompass,
  Thermometer,
} from 'svelte-widgets/icons'
import benchmark_counts from '$routes/data/benchmark-element-counts.json'
import wbm_counts from '$routes/data/wbm-element-counts-by-occurrence.json'

const n_elements = (counts: Record<string, number | null>) =>
  Object.values(counts).filter((count) => count !== null && count > 0).length

const { n_structures } = DATASETS.WBM
if (n_structures === null) throw new Error(`Missing WBM structure count: ${n_structures}`)
const file_url = (key: string) => {
  const file = data_files[key]
  if (!file || typeof file === `string` || !file.url)
    throw new Error(`No public download URL for ${key}`)
  return file.url
}

export const test_sets = {
  wbm: {
    id: `wbm`,
    name: `WBM`,
    credit: [`Wang, Botti, and Marques`, `https://doi.org/10.1038/s41524-020-00481-6`],
    coverage: `${format_num(n_structures, `,`)} crystal structures · ${n_elements(wbm_counts)} elements`,
    element_counts: wbm_counts,
    count_unit: `Structures`,
    description: `Chemically substituted crystals with DFT relaxations test stability predictions and how closely relaxed geometries reproduce DFT. Discovery uses a unique-prototype subset; geometry optimization evaluates the relaxed structures.`,
    availability: `Public initial and relaxed structures, energies, and symmetry references.`,
    links: [
      [`Dataset details`, `/data/wbm`, Info],
      [`Initial structures`, file_url(`wbm_initial_atoms`), Download],
      [`DFT-relaxed structures`, file_url(`wbm_relaxed_atoms`), Download],
    ],
  },
  phonondb: {
    id: `phonondb`,
    name: `PhononDB PBE`,
    credit: [
      `Togo`,
      `https://github.com/atztogo/phonondb#url-links-to-phono3py-finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe`,
    ],
    coverage: `103 crystal structures · ${n_elements(benchmark_counts.phonondb)} elements`,
    element_counts: benchmark_counts.phonondb,
    count_unit: `Structures`,
    description: `Materials Project crystals with PBE reference calculations test harmonic phonons and anharmonic lattice thermal conductivity at 300 K.`,
    availability: `Public structures and thermal-conductivity references without non-analytical corrections.`,
    links: [
      [`Structures`, file_url(`phonondb_pbe_103_structures`), Download],
      [`DFT references`, file_url(`phonondb_pbe_103_kappa_no_nac`), Download],
    ],
  },
  dynamat: {
    id: `dynamat`,
    name: `DynaMat v1.0`,
    credit: [`Gawkowski et al.`, `https://arxiv.org/abs/2607.03433`],
    coverage: `17 systems · ${n_elements(benchmark_counts.dynamat)} elements · 293–1500 K`,
    element_counts: benchmark_counts.dynamat,
    count_unit: `Systems`,
    description: `Ab-initio NVT trajectories span metals, alloys, perovskites, molecular crystals, and transition-metal dichalcogenides. They test structural, thermodynamic, and vibrational observables at finite temperature.`,
    availability: `Public trajectories omit energies and forces. These withheld labels are used only for maintainer-computed diagnostics, outside CMDS.`,
    links: [
      [`Trajectories (HDF5)`, file_url(`dynamat_v1_0_md_trajectories`), Download],
      [`Reference paper`, `https://arxiv.org/abs/2607.03433`, Info],
    ],
  },
  diatomics: {
    id: `diatomics`,
    name: `Diatomic curves`,
    credit: null,
    coverage: `92 homonuclear dimers · ${n_elements(benchmark_counts.diatomics)} elements · PBE and r2SCAN`,
    element_counts: benchmark_counts.diatomics,
    count_unit: `Dimers`,
    description: `Potential-energy and force curves from H₂ to U₂ probe bond geometry, dissociation, repulsion, and smoothness. PBE reference-relative errors exclude low-quality reference curves.`,
    availability: `Public DFT curves include energies, forces, and site-projected magnetic moments. Reference-relative scores use PBE.`,
    links: [
      [`DFT curves (.json.gz)`, file_url(`diatomics_dft_reference`), Download],
      [
        `Reference provenance`,
        `${pkg.repository}/blob/main/site/src/lib/diatomics-dft.readme.md`,
        Info,
      ],
    ],
  },
} as const
export type BenchmarkDataset = (typeof test_sets)[keyof typeof test_sets]

export const benchmarks = {
  discovery: { dataset: test_sets.wbm, icon: Search, key: `discovery` },
  'geo-opt': { dataset: test_sets.wbm, icon: RulerSquareCompass, key: `geo_opt` },
  phonons: { dataset: test_sets.phonondb, icon: Phonons, key: `phonons` },
  md: { dataset: test_sets.dynamat, icon: Thermometer, key: `md` },
  diatomics: { dataset: test_sets.diatomics, icon: Molecule, key: `diatomics` },
} as const
export type BenchmarkTask = keyof typeof benchmarks
