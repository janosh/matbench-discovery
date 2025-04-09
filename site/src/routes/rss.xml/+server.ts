import { MODELS } from '$lib'
import { format_train_set } from '$lib/metric-helpers'
import type { ModelData } from '$lib/types'
import { description, homepage, name } from '$site/package.json'
import { pretty_num } from 'elementari'

export const prerender = true

// Ensure homepage URL ends with a slash for consistent path joining
const base_url = homepage.endsWith(`/`) ? homepage : `${homepage}/`

// Formats model data as XML for RSS feed
function format_model_for_rss(model: ModelData): string {
  // Extract metrics from the 'full_test_set' or fallback to first available discovery set
  const discovery_metrics =
    model.metrics?.discovery?.full_test_set ||
    (model.metrics?.discovery && Object.values(model.metrics.discovery)[0])

  const training_set = format_train_set(model.training_set)
  // Remove HTML tags for plain text display
  const clean_training_set = training_set.replace(/<[^>]*>/g, ``)

  const metrics_text = discovery_metrics
    ? Object.entries(discovery_metrics)
        .filter(([_key, value]) => typeof value === `number`)
        .map(([key, value]) => `${key}: ${pretty_num(value as number)}`)
        .join(`, `)
    : `No metrics available`

  const authors_text = model.authors
    ? model.authors
        .map((author) => {
          const parts = []
          parts.push(author.name)
          if (author.affiliation) parts.push(`(${author.affiliation})`)
          return parts.join(` `)
        })
        .join(`, `)
    : `Unknown authors`

  return `
    <h2>Model: ${model.model_name}</h2>
    <p><strong>Authors:</strong> ${authors_text}</p>
    <p><strong>Date Added:</strong> ${model.date_added}</p>
    <p><strong>Training Set:</strong> ${clean_training_set}</p>
    <p><strong>Parameters:</strong> ${pretty_num(model.model_params)}</p>
    <p><strong>Metrics:</strong> ${metrics_text}</p>
    <p><strong>Targets:</strong> ${model.targets}</p>
    <p>
      <a href="${base_url}models/${model.model_key}">View model details</a>
      ${model.paper ? `| <a href="${model.paper}">Read paper</a>` : ``}
      ${model.repo ? `| <a href="${model.repo}">View code repository</a>` : ``}
    </p>
  `.trim()
}

/**
 * Generates an RSS feed of all models
 */
export async function GET() {
  // Sort models by date added (newest first)
  const sorted_models = [...MODELS].sort((a, b) => {
    return new Date(b.date_added).getTime() - new Date(a.date_added).getTime()
  })

  const headers = { 'Content-Type': `application/xml` }

  const xml = `
    <rss xmlns:atom="http://www.w3.org/2005/Atom" version="2.0">
      <channel>
        <title>${name} - Model Leaderboard</title>
        <description>${description}</description>
        <link>${base_url}</link>
        <atom:link href="${base_url}rss.xml" rel="self" type="application/rss+xml"/>
        ${sorted_models
          .map(
            (model) => `
            <item>
              <title>${model.model_name}</title>
              <description><![CDATA[${format_model_for_rss(model)}]]></description>
              <link>${base_url}models/${model.model_key}</link>
              <guid isPermaLink="true">${base_url}models/${model.model_key}</guid>
              <pubDate>${new Date(model.date_added).toUTCString()}</pubDate>
            </item>
          `,
          )
          .join(``)}
      </channel>
    </rss>
  `.trim()

  return new Response(xml, { headers })
}
