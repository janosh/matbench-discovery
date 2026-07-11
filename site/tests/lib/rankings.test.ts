import { ALL_METRICS } from '$lib/labels'
import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
import type { Label, ModelData } from '$lib/types'
import { describe, expect, it } from 'vitest'

// minimal model stub with just the fields the ranked metric paths read
const make_model = (
  model_key: string,
  { cps, f1 }: { cps?: number; f1?: number } = {},
): ModelData =>
  ({
    model_key,
    CPS: cps ?? Number.NaN,
    metrics: { discovery: { unique_prototypes: { F1: f1 } } },
  }) as unknown as ModelData

const models = [
  make_model(`model-a`, { cps: 0.9, f1: 0.8 }),
  make_model(`model-b`, { cps: 0.7, f1: 0.9 }),
  make_model(`model-c`, { cps: 0.8, f1: 0.8 }),
]

const rank_of = (model_key: string, metric_key: string) =>
  model_metric_ranks(model_key, models, RANKED_METRICS).find(
    ({ metric }) => metric.key === metric_key,
  )
const rank_keys = (model_key: string, metrics: readonly Label[] = RANKED_METRICS) =>
  model_metric_ranks(model_key, models, metrics).map(({ metric }) => metric.key)

describe(`model_metric_ranks`, () => {
  it.each([
    [`model-a`, `CPS`, 1, 3], // higher=better
    [`model-b`, `CPS`, 3, 3],
  ])(`%s ranks %s #%d of %d`, (model_key, metric_key, rank, n_models) => {
    expect(rank_of(model_key, metric_key)).toMatchObject({ rank, n_models })
  })

  it(`ties share the better rank`, () => {
    expect(rank_of(`model-a`, `F1`)).toMatchObject({ rank: 2, n_models: 3 })
    expect(rank_of(`model-c`, `F1`)).toMatchObject({ rank: 2, n_models: 3 })
    expect(rank_of(`model-b`, `F1`)).toMatchObject({ rank: 1, n_models: 3 })
  })

  it(`omits metrics the model has no value for`, () => {
    expect(rank_keys(`model-c`)).toEqual([`CPS`, `F1`])
  })

  it(`returns [] for unknown model keys`, () => {
    expect(model_metric_ranks(`no-such-model`, models, RANKED_METRICS)).toEqual([])
  })

  it(`custom metric list restricts which metrics get ranked`, () => {
    expect(rank_keys(`model-a`, [ALL_METRICS.F1])).toEqual([`F1`])
  })

  it(`RANKED_METRICS carry task-prefixed labels and leaderboard hrefs`, () => {
    expect(
      RANKED_METRICS.map(({ label, rank_href }) => [label, rank_href]),
    ).toStrictEqual([
      [`CPS`, `/`],
      [`Discovery F1`, `/tasks/discovery`],
      [`Geo Opt RMSD`, `/tasks/geo-opt`],
      [`Phonons κ<sub>SRME</sub>`, `/tasks/phonons`],
      [`MD CMDS`, `/tasks/md`],
      [`Diatomics CDS`, `/tasks/diatomics`],
    ])
  })
})
