// This file is auto-generated from dataset-schema.yml. Do not edit directly.

/**
 * Computational or experimental method used to generate the data
 */
export type MethodEnum = `DFT` | `experiment` | `ML`
/**
 * DFT code used for calculations
 */
export type DftCodeEnum =
  | `VASP`
  | `Quantum ESPRESSO`
  | `CASTEP`
  | `AbInit`
  | `FHI-aims`
  | `CP2K`
  | `Various`
/**
 * Exchange-correlation functional used in DFT calculations
 */
export type FunctionalEnum =
  | `PBE`
  | `PBEsol`
  | `PBE+U`
  | `SCAN`
  | `SCAN+U`
  | `r2SCAN`
  | `r2SCAN+U`
  | `HSE`
  | `HSE06`
  | `PBE0`
  | `Various`
/**
 * Pseudopotentials used in DFT calculations
 */
export type PseudoPotentialsEnum = `PBE` | `PBE_52` | `PBE_54` | `PBE_64` | `Various`

export interface DatasetRecord {
  [k: string]: Dataset
}
/**
 * This interface was referenced by `DatasetRecord`'s JSON-Schema definition
 * via the `patternProperty` "^[a-zA-Z0-9 ]+$".
 */
export interface Dataset {
  /**
   * Full name of the dataset
   */
  title: string
  /**
   * Version of the dataset
   */
  version?: string
  /**
   * Detailed description of the dataset
   */
  description: string
  /**
   * Number of structures in the dataset
   */
  n_structures: number
  /**
   * Number of unique materials in the dataset
   */
  n_materials?: number
  /**
   * Chemical elements included in the dataset
   */
  elements?: string[]
  /**
   * Temperature range of structures in the dataset (e.g., '0-5000 K')
   */
  temperature_range?: string
  /**
   * Pressure range of structures in the dataset (e.g., '0-1000 GPa')
   */
  pressure_range?: string
  /**
   * Primary URL for the dataset
   */
  url: string
  /**
   * URL to download the dataset
   */
  download_url?: string
  /**
   * DOI reference for the dataset
   */
  doi: string
  /**
   * Other datasets this dataset is derived from
   */
  derived_from?: string[]
  /**
   * Whether the dataset is openly available
   */
  open: boolean
  /**
   * License under which the dataset is published
   */
  license:
    | `MIT`
    | `Apache-2.0`
    | `CC-BY-4.0`
    | `GPL-3.0`
    | `BSD-3-Clause`
    | `LGPL-3.0`
    | `Meta Research`
    | `unreleased`
  /**
   * People or organizations who created the dataset
   */
  created_by?: Person[]
  /**
   * Date when the dataset was created
   */
  date_created: string
  /**
   * Date when the dataset was added to this collection
   */
  date_added?: string
  /**
   * Parameters and methods used to generate the dataset
   */
  params?: {
    /**
     * Method(s) used to generate the data
     */
    method?: MethodEnum | MethodEnum[]
    /**
     * DFT code(s) used for calculations
     */
    code?: DftCodeEnum | DftCodeEnum[]
    /**
     * Version of the DFT code used
     */
    code_version?: string
    /**
     * Exchange-correlation functional(s) used
     */
    functional?: FunctionalEnum | FunctionalEnum[]
    /**
     * Pseudopotential(s) used
     */
    pseudopotentials?: PseudoPotentialsEnum | PseudoPotentialsEnum[]
    [k: string]: unknown
  }
}
export interface Person {
  name: string
  affiliation?: string
  email?: string
  url?: string
  orcid?: string
  github?: string
  corresponding?: boolean
}
