<script lang="ts">
  import { goto } from '$app/navigation'
  import { spacegroup_html, type KappaMaterialModelRow } from '#lib/materials.js'
  import {
    as_phonon_dos,
    has_kappa_parity_model,
    load_kappa_parity_base,
    load_kappa_parity_model,
  } from '#lib/parity/kappa-parity.js'
  import {
    build_phonon_mode_data,
    has_kappa_modes,
    load_kappa_modes,
  } from '#lib/parity/kappa-modes.js'
  import { labeled_vline } from '#lib/fig-helpers.js'
  import ModelSelect from '#lib/ModelSelect.svelte'
  import { format_num } from 'matterviz/labels'
  import { ScatterPlot, type DataSeries } from 'matterviz/plot'
  import { sanitize_compact_formula } from 'matterviz/sanitize'
  import type { PhononDos } from 'matterviz/spectral'
  import { spacegroup_to_crystal_sys } from 'matterviz/symmetry'
  import { HeatmapTable, type Column } from 'matterviz/table'
  import { MultiSelect, Spinner } from 'svelte-widgets'
  import { to_error } from 'svelte-widgets/utils'
  import type { PageData } from './$types'

  let { data }: { data: PageData } = $props()
  let material = $derived(data.material)
  let { material_id, formula, spg_num, kappa_dft } = $derived(material)
  let formula_html = $derived(sanitize_compact_formula(formula))
  let spg_label = $derived(
    `${spacegroup_html(spg_num)}, ${spacegroup_to_crystal_sys(spg_num) ?? `unknown system`}`,
  )

  const failure_labels = [
    [`imag_modes`, `imaginary modes`],
    [`broken_sym`, `broken symmetry`],
    [`max_steps`, `relaxation hit max steps`],
    [`below_floor`, `κ below half the diffuson minimum`],
  ] as const
  const failures = (row: KappaMaterialModelRow): string =>
    failure_labels.flatMap(([key, label]) => (row[key] ? [label] : [])).join(`, `) ||
    (row.srme_censored ? `κ calculation failed` : ``)

  let scored_rows = $derived(
    material.rows.filter(
      (row): row is KappaMaterialModelRow & { kappa_ml: number; srme: number } =>
        row.kappa_ml !== null && row.kappa_ml > 0 && row.srme !== null,
    ),
  )
  // models whose κ lies within a factor 2, ±30% and ±10% of the DFT value
  const within_dft = [
    (ratio: number) => ratio >= 0.5 && ratio <= 2,
    (ratio: number) => Math.abs(ratio - 1) <= 0.3,
    (ratio: number) => Math.abs(ratio - 1) <= 0.1,
  ]
  let n_within = $derived(
    within_dft.map(
      (within) =>
        material.rows.filter(
          ({ kappa_ratio }) => kappa_ratio !== null && within(kappa_ratio),
        ).length,
    ),
  )

  // === κ_ML vs κ_SRME scatter ===
  let series = $derived<DataSeries<KappaMaterialModelRow & Record<string, unknown>>[]>([
    {
      x: scored_rows.map((row) => row.kappa_ml),
      y: scored_rows.map((row) => row.srme),
      metadata: scored_rows,
      markers: `points`,
      // censored SRMEs (failed predictions) render hollow, like on the phonons page
      point_style: scored_rows.map((row) => ({
        fill: row.srme_censored ? `transparent` : `#636efa`,
        stroke: row.srme_censored ? `#ef553b` : `white`,
        stroke_width: row.srme_censored ? 1.5 : 0.5,
        radius: 6,
      })),
    },
  ])
  let kappa_range = $derived.by((): [number, number] => {
    const values = [...scored_rows.map((row) => row.kappa_ml), kappa_dft ?? 1]
    return [Math.min(...values) / 1.5, Math.max(...values) * 1.5]
  })

  // === per-model table ===
  let table_rows = $derived(
    material.rows.map((row) => ({
      Model: `<a href="/models/${row.model_key}">${row.label}</a>`,
      kappa_ml: row.kappa_ml,
      kappa_ratio: row.kappa_ratio,
      srme: row.srme,
      freq_w1: row.freq_w1,
      Issues: failures(row),
    })),
  )
  const all_columns: Column<(typeof table_rows)[number]>[] = [
    { id: `Model`, label: `Model`, pin: `left` },
    {
      id: `srme`,
      label: `κ<sub>SRME</sub>`,
      format: `.3f`,
      better: `lower`,
      description: `Per-material mode-resolved κ error (0 = perfect, 2 = maximum error or failure)`,
    },
    {
      id: `kappa_ml`,
      label: `κ<sub>ML</sub>`,
      format: `.4~`,
      description: `Predicted lattice thermal conductivity at 300 K (W/mK)`,
    },
    {
      id: `kappa_ratio`,
      label: `κ<sub>ML</sub>/κ<sub>DFT</sub>`,
      format: `.3~`,
      description: `1 = exact; below 1 = underestimated`,
    },
    {
      id: `freq_w1`,
      label: `W<sub>1</sub>(ω)`,
      format: `.3f`,
      better: `lower`,
      description: `Wasserstein-1 distance between ML and DFT phonon frequency spectra (THz)`,
    },
    { id: `Issues`, label: `Issues` },
  ]
  // the Issues column only when some model failed on this material
  let columns = $derived(
    table_rows.some(({ Issues }) => Issues)
      ? all_columns
      : all_columns.filter(({ id }) => id !== `Issues`),
  )

  // === phonon DOS overlay (assets load on demand: ~60 kB per model) ===
  const model_option = ({ label, model_key }: KappaMaterialModelRow) => ({
    label,
    value: model_key,
  })
  let dos_options = $derived(
    material.rows
      .filter((row) => has_kappa_parity_model(row.model_key))
      .map(model_option),
  )
  // default: this material's best and worst scored model, so the spread shows at once
  let dos_selected = $derived.by(() => {
    const keys = new Set([material.best?.model_key, material.worst?.model_key])
    return dos_options.filter(({ value }) => keys.has(value))
  })
  const load_spectral = () => import(`matterviz/spectral`)
  // assets are client fetches: prerender must not start them (SSR evaluates {#await})
  const browser = !import.meta.env.SSR
  const load_doses = async (keys: string[]) => {
    const [spectral, base, ...models] = await Promise.all([
      load_spectral(),
      load_kappa_parity_base(),
      ...keys.map(load_kappa_parity_model),
    ])
    const doses: Record<string, PhononDos> = {}
    const dft = as_phonon_dos(base.dft_dos[material_id])
    if (dft) doses[`DFT (PBE)`] = dft
    for (const [idx, model] of models.entries()) {
      const dos = as_phonon_dos(model.ml_dos[material_id])
      if (dos) doses[dos_selected[idx]?.label ?? model.model_key] = dos
    }
    return { spectral, doses }
  }
  let doses_promise = $derived(
    browser ? load_doses(dos_selected.map((option) => String(option.value))) : null,
  )

  // === phonon band structure + modes (~1 MB per model; defaults to the best model) ===
  let modes_options = $derived(
    material.rows.filter((row) => has_kappa_modes(row.model_key)).map(model_option),
  )
  let modes_selected = $derived(modes_options[0] ?? null)
  let modes_key = $derived(modes_selected && String(modes_selected.value))
  let modes_promise = $derived(
    browser && modes_key && Promise.all([load_spectral(), load_kappa_modes(modes_key)]),
  )
</script>

<svelte:head>
  <title>{material_id} {formula} • Phonon materials • Matbench Discovery</title>
</svelte:head>

<h1>{@html formula_html} <small>{material_id}</small></h1>

<dl class="facts">
  <div>
    <dt>Spacegroup</dt>
    <dd>{@html spg_label}</dd>
  </div>
  <div>
    <dt>PBE κ<sub>L</sub> (300 K)</dt>
    <dd>{kappa_dft === null ? `n/a` : `${format_num(kappa_dft, `.4~`)} W/mK`}</dd>
  </div>
  <div>
    <dt>Median κ<sub>SRME</sub></dt>
    <dd>
      {material.median_srme === null ? `n/a` : format_num(material.median_srme, `.3f`)}
      {#if material.median_srme_rank !== null}
        <small title="Rank among all materials, from lowest (1) to highest median error"
          >#{material.median_srme_rank}/{material.n_ranked}</small
        >
      {/if}
    </dd>
  </div>
  <div>
    <dt>Within 2× · 30% · 10% of DFT κ</dt>
    <dd>
      {n_within.join(` · `)}
      <small>of {material.n_models} models</small>
    </dd>
  </div>
  <div>
    <dt>Failed predictions</dt>
    <dd>{material.n_censored} of {material.n_models} models</dd>
  </div>
</dl>

<section>
  <h2>Predicted κ vs per-material error</h2>
  <p class="caption">
    Each point is one model. The dashed line marks the PBE reference.
    {#if scored_rows.some(({ srme_censored }) => srme_censored)}Hollow red markers are
      failed predictions censored to κ<sub>SRME</sub> = 2.{/if}
    κ<sub>SRME</sub> sums errors over individual phonon modes, so it can only exceed the error
    of the total κ: a point near the dashed line but high up gets the total right through error
    cancellation, with some modes conducting too much heat and others too little. Click a point
    to open the model.
  </p>
  <ScatterPlot
    {series}
    style="height: 420px"
    ref_lines={kappa_dft ? [labeled_vline(kappa_dft, `PBE κ`)] : []}
    x_axis={{
      label: `Predicted κ<sub>L</sub> (W/mK)`,
      scale_type: `log`,
      format: `~s`,
      range: kappa_range,
    }}
    y_axis={{ label: `κ<sub>SRME</sub>`, range: [0, 2.05], format: `.1f` }}
    on_point_click={({ point }) => {
      const row = scored_rows[point.point_idx]
      if (row) void goto(`/models/${row.model_key}`)
    }}
  >
    {#snippet tooltip({ metadata })}
      {#if metadata}
        {@const row = metadata as KappaMaterialModelRow}
        <strong>{row.label}</strong><br />
        κ: {format_num(row.kappa_ml ?? NaN, `.3~`)} <small>W/mK</small>
        {#if row.kappa_ratio !== null}(×{format_num(row.kappa_ratio, `.2~`)} DFT){/if}<br
        />
        κ<sub>SRME</sub>: {format_num(row.srme ?? NaN, `.3~`)}
        {#if failures(row)}<br /><em>{failures(row)}</em>{/if}
      {/if}
    {/snippet}
  </ScatterPlot>
</section>

<section>
  <h2>All models</h2>
  <HeatmapTable
    data={table_rows}
    {columns}
    sort={[{ column: `srme`, dir: `asc` }]}
    default_num_format=".3~f"
  />
</section>

<section>
  <h2>Phonon density of states</h2>
  <p class="caption">
    Harmonic phonon DOS of the DFT reference and the selected models, computed from each
    model's force constants on its own relaxed structure. Defaults to this material's best
    and worst model.
  </p>
  <div class="dos-models" style="display: flex; justify-content: center">
    <ModelSelect
      options={dos_options}
      bind:value={dos_selected}
      placeholder="Add models to compare"
    />
  </div>
  {#if doses_promise}
    {#await doses_promise}
      <Spinner text="Loading phonon DOS..." />
    {:then { spectral, doses }}
      {#if Object.keys(doses).length}
        <spectral.Dos {doses} style="height: 400px" padding={{ t: 20, b: 60, r: 10 }} />
      {:else}
        <p class="caption">No phonon DOS stored for this material.</p>
      {/if}
    {:catch error}
      <p class="caption" role="alert">
        Could not load phonon data: {to_error(error).message}
      </p>
    {/await}
  {/if}
</section>

{#if modes_options.length}
  <section>
    <h2>Phonon band structure and modes</h2>
    <div class="modes-controls">
      <MultiSelect
        options={modes_options}
        mode="single"
        min_select={1}
        bind:value={modes_selected}
        style="width: 24em; max-width: 100%; border: 1px solid var(--border)"
      />
    </div>
    {#if modes_promise}
      {#await modes_promise}
        <Spinner text="Loading phonon modes..." />
      {:then [spectral, modes_model]}
        {@const modes = modes_model.materials[material_id]}
        {#if modes}
          <spectral.PhononModeExplorer
            dataset={{
              modes: build_phonon_mode_data(modes),
              filename: `${material_id}-${modes_model.model_key}`,
            }}
            style="height: 520px"
            supercell={[2, 2, 2]}
          />
        {:else}
          <p class="caption">No phonon modes stored for this material.</p>
        {/if}
      {:catch error}
        <p class="caption" role="alert">
          Could not load phonon modes: {to_error(error).message}
        </p>
      {/await}
    {/if}
  </section>
{/if}

<nav class="material-nav" aria-label="Material navigation">
  {#if material.prev_id}<a href="/materials/phonons/{material.prev_id}"
      >&larr; {material.prev_id}</a
    >{:else}<span></span>{/if}
  <a href="/materials/phonons">All materials</a>
  {#if material.next_id}<a href="/materials/phonons/{material.next_id}"
      >{material.next_id} &rarr;</a
    >{:else}<span></span>{/if}
</nav>

<style>
  h1 small {
    font-size: 0.5em;
    color: var(--text-secondary);
  }
  .material-nav {
    display: flex;
    justify-content: space-between;
    margin-block: 2.5em 1em;
  }
  .facts {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 0.8em 2.5em;
    div {
      text-align: center;
    }
    dt {
      font-size: 0.85em;
      color: var(--text-secondary);
    }
    dd {
      margin: 0;
      font-weight: 600;
      small {
        font-weight: normal;
        color: var(--text-secondary);
      }
    }
  }
  section {
    margin-block: 2.5em;
  }
  h2 {
    text-align: center;
  }
  .caption {
    text-align: center;
    color: var(--text-secondary);
    max-width: 50em;
    margin-inline: auto;
  }
  .modes-controls {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 0.8em;
    margin-block: 1em;
  }
</style>
