#!/usr/bin/env -S deno run -A
// deno-lint-ignore-file no-await-in-loop
// Fetch GitHub activity (stars, forks, commits, contributors) for MLIP repos from models/*.yml
// Usage: deno run -A scripts/fetch-mlip-repos-github-activity.ts [--force-refresh]
// Set GITHUB_TOKEN env var to avoid rate limits
import { parse as parseYAML } from 'jsr:@std/yaml@^1.0.5'
import { mkdir, readdir, readFile, stat, writeFile } from 'node:fs/promises'
import { dirname, join } from 'node:path'
import process from 'node:process'
import { fileURLToPath } from 'node:url'

const __dirname = dirname(fileURLToPath(import.meta.url))
const site_root = join(__dirname, `..`)
const project_root = join(site_root, `..`)
const cache_dir = join(site_root, `.cache`)
const output_file = join(site_root, `src/routes/models/mlip-github-activity.json`)
const models_dir = join(project_root, `models`)

const CACHE_TTL_MS = 24 * 60 * 60 * 1000 // 24 hours
const force_refresh = process.argv.includes(`--force-refresh`)

type RepoData = {
  name: string
  repo: string
  stars: number
  forks: number
  commits_last_year: number
  contributors: number
}

const extract_github_repo = (url: string | null): string | null => {
  if (!url) return null
  const match = /github\.com\/([^/]+\/[^/]+)/.exec(url)
  return match ? match[1].replace(/\.git$/, ``) : null
}

const load_repos_from_models = async (): Promise<{ name: string; repo: string }[]> => {
  const repos_map = new Map<string, string>() // repo -> name

  const model_dirs = await readdir(models_dir, { withFileTypes: true })

  for (const dir_entry of model_dirs) {
    if (!dir_entry.isDirectory()) continue

    const model_dir = join(models_dir, dir_entry.name)
    const files = await readdir(model_dir)
    const yml_files = files.filter((file) => file.endsWith(`.yml`))

    for (const yml_file of yml_files) {
      try {
        const yml_content = await readFile(join(model_dir, yml_file), `utf-8`)
        const data = parseYAML(yml_content) as { model_name?: string; repo?: string }

        const repo_url = extract_github_repo(data.repo ?? null)
        if (repo_url && data.model_name && !repos_map.has(repo_url)) {
          const label = data.model_name.split(`-`)[0] || data.model_name
          repos_map.set(repo_url, label)
        }
      } catch (err) {
        console.warn(`  Warning: ${yml_file} - ${(err as Error).message}`)
      }
    }
  }

  return Array.from(repos_map.entries()).map(([repo, name]) => ({ name, repo }))
}

const get_cached = async (cache_file: string): Promise<RepoData | null> => {
  if (force_refresh) return null
  try {
    const age_ms = Date.now() - (await stat(cache_file)).mtimeMs
    if (age_ms < CACHE_TTL_MS) {
      const cached = JSON.parse(await readFile(cache_file, `utf-8`))
      if (cached.name && cached.repo && cached.stars !== undefined) {
        return cached as RepoData
      }
    }
  } catch {
    // Cache miss or invalid format
  }
  return null
}

const fetch_github = async (url: string, headers: Record<string, string>) => {
  const res = await fetch(url, { headers })
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  return res
}

const get_count_from_pagination = (link: string | null) => {
  const match = /page=(\d+)>; rel="last"/.exec(link || ``)
  return match ? parseInt(match[1]) : null
}

const get_github_stats = async (
  name: string,
  repo: string,
  token: string | null,
): Promise<RepoData | null> => {
  const cache_file = join(cache_dir, `${repo.replace(`/`, `_`)}.json`)
  const cached = await get_cached(cache_file)
  if (cached) {
    console.log(`✓ ${repo} (cached)`)
    return cached
  }

  const headers: Record<string, string> = { Accept: `application/vnd.github+json` }
  if (token) headers.Authorization = `Bearer ${token}`

  try {
    // Fetch repo info
    const repo_res = await fetch_github(`https://api.github.com/repos/${repo}`, headers)
    const repo_data = await repo_res.json()

    // Fetch commit count (last year)
    const year_ago = new Date()
    year_ago.setFullYear(year_ago.getFullYear() - 1)
    const commits_res = await fetch_github(
      `https://api.github.com/repos/${repo}/commits?since=${year_ago.toISOString()}&per_page=1`,
      headers,
    )
    const commits_last_year =
      get_count_from_pagination(commits_res.headers.get(`Link`)) ??
        ((await commits_res.json()).length || 0)

    // Fetch contributor count
    const contrib_res = await fetch_github(
      `https://api.github.com/repos/${repo}/contributors?per_page=1&anon=true`,
      headers,
    )
    const contributors = get_count_from_pagination(contrib_res.headers.get(`Link`)) ??
      ((await contrib_res.json()).length || 0)
    const stars = repo_data.stargazers_count || 0
    const forks = repo_data.forks_count || 0
    const data: RepoData = { name, repo, stars, forks, commits_last_year, contributors }

    // Cache the result
    await mkdir(cache_dir, { recursive: true })
    await writeFile(cache_file, JSON.stringify(data))
    console.log(
      `✓ ${repo}: ${data.stars}★ ${data.forks}f ${data.commits_last_year}c ${data.contributors}contrib`,
    )
    return data
  } catch (err) {
    console.error(`✗ ${repo}:`, (err as Error).message)
    return null
  }
}

const main = async () => {
  let existing_map = new Map<string, RepoData>()
  try {
    const data = JSON.parse(await readFile(output_file, `utf-8`)) as RepoData[]
    existing_map = new Map(data.map((entry) => [entry.repo, entry]))
  } catch {
    // File doesn't exist yet, start fresh
  }

  const repos = await load_repos_from_models()
  const token = process.env.GITHUB_TOKEN ?? null
  const results: RepoData[] = []

  for (const { name, repo } of repos) {
    const existing = existing_map.get(repo)
    if (existing && !force_refresh) {
      const cache_file = join(cache_dir, `${repo.replace(`/`, `_`)}.json`)
      try {
        const age_ms = Date.now() - (await stat(cache_file)).mtimeMs
        if (age_ms < CACHE_TTL_MS) {
          results.push(existing)
          continue
        }
      } catch {
        // Cache file missing, refetch from GitHub
      }
    }

    const stats = await get_github_stats(name, repo, token)
    if (stats) results.push(stats)
    else if (existing) results.push(existing)
    await new Promise((resolve) => setTimeout(resolve, 300))
  }

  await writeFile(output_file, JSON.stringify(results, null, 2))
  console.log(`Saved ${results.length}/${repos.length} repos`)
}

main().catch((err) => {
  console.error(`Fatal:`, err)
  writeFile(output_file, JSON.stringify([])).catch(() => {})
})
