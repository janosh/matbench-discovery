#!/usr/bin/env node
// Fetch GitHub activity data for MLIP repositories during site builds

// deno-lint-ignore-file no-await-in-loop
import { mkdir, readFile, writeFile } from 'fs/promises'
import process from 'node:process'
import { dirname, join } from 'path'
import { fileURLToPath } from 'url'

const __dirname = dirname(fileURLToPath(import.meta.url))
const site_root = join(__dirname, `..`)
const cache_dir = join(site_root, `.cache`)
const output_file = join(site_root, `src/lib/github-activity.json`)
const mlip_models_path = join(
  site_root,
  `../../papers/foundation-models-chem-rev/figs/mlip_models.json`,
)

// Extract last page number from GitHub API Link header
function get_last_page(link_header) {
  const match = /page=(\d+)>; rel="last"/.exec(link_header || ``)
  return match ? parseInt(match[1]) : null
}

// Get GitHub stats with caching
async function get_github_stats(repo, token = null) {
  const cache_file = join(cache_dir, `${repo.replace(`/`, `_`)}.json`)

  // Check cache first
  try {
    const cached = await readFile(cache_file, `utf-8`)
    console.log(`✓ Using cached data for ${repo}`)
    return JSON.parse(cached)
  } catch {
    // Cache miss, fetch from API
  }

  const headers = { Accept: `application/vnd.github.v3+json` }
  if (token) headers.Authorization = `Bearer ${token}`

  try {
    const response = await fetch(`https://api.github.com/repos/${repo}`, { headers })
    if (!response.ok) {
      console.error(`Error fetching ${repo}: ${response.statusText}`)
      return null
    }

    const data = await response.json()
    const one_year_ago = new Date()
    one_year_ago.setFullYear(one_year_ago.getFullYear() - 1)

    // Get commits in the last year
    const commits_res = await fetch(
      `https://api.github.com/repos/${repo}/commits?since=${one_year_ago.toISOString()}&per_page=100`,
      { headers },
    )
    if (commits_res.ok) {
      const commits = await commits_res.json()
      const last_page = get_last_page(commits_res.headers.get(`Link`))
      data.commits_last_year = last_page
        ? (last_page - 1) * 100 + commits.length
        : commits.length
    } else {
      console.warn(`Warning: Could not fetch commits for ${repo}`)
      data.commits_last_year = 0
    }

    // Get number of contributors
    const contrib_res = await fetch(
      `https://api.github.com/repos/${repo}/contributors?per_page=1&anon=true`,
      { headers },
    )
    if (contrib_res.ok) {
      const last_page = get_last_page(contrib_res.headers.get(`Link`))
      data.contributors = last_page || (await contrib_res.json()).length
    } else {
      console.warn(`Warning: Could not fetch contributors for ${repo}`)
      data.contributors = 0
    }

    // Cache the result
    await mkdir(cache_dir, { recursive: true })
    await writeFile(cache_file, JSON.stringify(data)).catch(() => {})

    return data
  } catch (error) {
    console.error(`Error fetching ${repo}:`, error.message)
    return null
  }
}

// Main function
async function main() {
  console.log(`Fetching GitHub repository data...`)

  // Read mlip_models.json
  let mlip_models
  try {
    mlip_models = JSON.parse(await readFile(mlip_models_path, `utf-8`))
  } catch (error) {
    console.error(
      `Error reading mlip_models.json from ${mlip_models_path}:`,
      error.message,
    )
    console.log(`Creating empty github-activity.json`)
    await writeFile(output_file, JSON.stringify([]))
    return
  }

  // Extract repositories
  const repos = Object.values(mlip_models)
    .filter(([, repo]) => repo)
    .map(([name, repo]) => ({ name, repo }))

  const token = process.env.GITHUB_TOKEN
  const repo_data = []

  // TODO might need to add delay here to avoid rate limiting
  for (const { name, repo } of repos) {
    const stats = await get_github_stats(repo, token)
    if (stats) {
      console.log(
        `✓ ${repo}: ${stats.stargazers_count || 0} stars, ${
          stats.forks_count || 0
        } forks, ${stats.commits_last_year || 0} commits, ${
          stats.contributors || 0
        } contributors`,
      )
      repo_data.push({
        name,
        repo,
        stars: stats.stargazers_count || 0,
        forks: stats.forks_count || 0,
        commits_last_year: stats.commits_last_year || 0,
        contributors: stats.contributors || 0,
      })
    }
  }

  await writeFile(output_file, JSON.stringify(repo_data, null, 2))
  console.log(`✓ GitHub activity data saved to ${output_file}`)
}

main().catch(console.error)
