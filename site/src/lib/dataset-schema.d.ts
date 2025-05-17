// This file is auto-generated from dataset-schema.yml. Do not edit directly.

/**
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "http_url".
 */
export type HttpUrl = string
/**
 * Computational or experimental method used to generate the data
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "method_enum".
 */
export type MethodEnum = `DFT` | `experiment` | `ML` | `GW` | `DMFT` | `MD`
/**
 * DFT code used for calculations
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "dft_code_enum".
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
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "functional_enum".
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
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "pseudo_potentials_enum".
 */
export type PseudoPotentialsEnum = `PBE` | `PBE_52` | `PBE_54` | `PBE_64` | `Various`
/**
 * License type:
 * - MIT: Massachusetts Institute of Technology
 * - CC-BY-4.0: Creative Commons Attribution 4.0 International
 * - CC-BY-NC-4.0: Creative Commons Attribution-NonCommercial 4.0 International
 * - CC-BY-SA-4.0: Creative Commons Attribution-ShareAlike 4.0 International
 * - GPL-3.0: GNU General Public v3.0
 * - BSD-3-Clause: Berkeley Software Distribution 3-Clause
 * - LGPL-3.0: GNU Lesser General Public License v3.0
 * - ASL: Academic Software License
 * - unreleased: No license since not released
 *
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "license_enum".
 */
export type LicenseEnum =
  | `MIT`
  | `Apache-2.0`
  | `CC-BY-4.0`
  | `CC-BY-SA-4.0`
  | `CC-BY-NC-4.0`
  | `GPL-3.0`
  | `BSD-3-Clause`
  | `LGPL-3.0`
  | `Meta Research`
  | `ASL`
  | `unreleased`

export interface DatasetRecord {
  [k: string]: Dataset
}
/**
 * This interface was referenced by `DatasetRecord`'s JSON-Schema definition
 * via the `patternProperty` "^[a-zA-Z0-9 ]+$".
 *
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "dataset".
 */
export interface Dataset {
  /**
   * Full name of the dataset
   */
  name: string
  /**
   * Slugified name of the dataset
   */
  slug?: string
  /**
   * Version of the dataset
   */
  version?: string
  /**
   * Detailed description of the dataset
   */
  description: string
  /**
   * HTML version of the dataset description
   */
  description_html?: string
  notes?: {
    [k: string]: unknown
  }
  /**
   * Number of structures in the dataset
   */
  n_structures: number
  /**
   * Number of unique materials in the dataset
   */
  n_materials?: number
  /**
   * Chemical elements included in the dataset (either symbols or atomic numbers)
   */
  elements?: (string | number)[]
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
   * Other datasets this dataset contains or is derived from
   */
  contains?: string[]
  /**
   * Whether the dataset is openly available
   */
  open: boolean
  /**
   * Whether the dataset is a static release or dynamically updated
   */
  static?: boolean
  /**
   * License under which the dataset is published
   */
  license:
    | `MIT`
    | `Apache-2.0`
    | `CC-BY-4.0`
    | `CC-BY-SA-4.0`
    | `CC-BY-NC-4.0`
    | `GPL-3.0`
    | `BSD-3-Clause`
    | `LGPL-3.0`
    | `Meta Research`
    | `ASL`
    | `unreleased`
  /**
   * URL to the OPTIMADE API endpoint or null if none
   */
  optimade_api?: null | HttpUrl
  /**
   * URL to the native API endpoint/documentation or null if none
   */
  native_api?: null | HttpUrl
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
   * Method(s) used to generate the data
   */
  method?: MethodEnum | MethodEnum[]
  /**
   * Parameters and methods used to generate the dataset
   */
  params?: {
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
/**
 * This interface was referenced by `DatasetRecord`'s JSON-Schema
 * via the `definition` "person".
 */
export interface Person {
  name: string
  affiliation?: string
  email?: string
  url?: string
  orcid?: string
  github?: string
  corresponding?: boolean
}
