// Auto-generates mlip-github-activity.json from models/*.yml (non-superseded only)
// Usage: npx tsx scripts/fetch-mlip-repos-github-activity.ts [--force-refresh]
// Set GITHUB_TOKEN env var to avoid rate limits
import { load as parseYAML } from 'js-yaml'
import { mkdir, readdir, readFile, stat, writeFile } from 'node:fs/promises'
import { join } from 'node:path'
import process from 'node:process'

const site_root = join(import.meta.dirname, `..`)
const project_root = join(site_root, `..`)
const cache_dir = join(site_root, `.cache`)
const output_file = join(site_root, `src/routes/models/mlip-github-activity.json`)
const models_dir = join(project_root, `models`)

const CACHE_TTL_MS = 24 * 60 * 60 * 1000 // 24 hours
const force_refresh = process.argv.includes(`--force-refresh`)

type RepoData = {
  name: string
  model_key: string
  repo: string
  stars: number
  forks: number
  commits_last_year: number
  contributors: number
}

const extract_github_repo = (url: string | null): string | null => {
  if (!url) return null
  const match = /github\.com\/(?<repo>[^/]+\/[^/]+)/.exec(url)
  return match?.groups?.repo.replace(/\.git$/, ``) ?? null
}

type ModelInfo = { name: string; model_key: string; repo: string }

const load_repos_from_models = async (): Promise<ModelInfo[]> => {
  // Map repo -> model info, preferring non-superseded models
  const repos_map = new Map<string, ModelInfo>()

  const model_dirs = (await readdir(models_dir, { withFileTypes: true }))
    .filter((d) => d.isDirectory())
    .sort((a, b) => a.name.localeCompare(b.name))

  for (const dir_entry of model_dirs) {
    const model_dir = join(models_dir, dir_entry.name)
    const files = await readdir(model_dir)
    const yml_files = files
      .filter((file) => file.endsWith(`.yml`))
      .sort((a, b) => a.localeCompare(b))

    for (const yml_file of yml_files) {
      try {
        const yml_content = await readFile(join(model_dir, yml_file), `utf-8`)
        const data = parseYAML(yml_content) as {
          model_name?: string
          model_key?: string
          repo?: string
          lifecycle?: string
        }

        if (data.lifecycle !== `active`) continue

        const repo_url = extract_github_repo(data.repo ?? null)
        if (repo_url && data.model_name && data.model_key && !repos_map.has(repo_url)) {
          // Only add if repo not already present (first non-superseded model wins)
          const label = data.model_name.split(`-`)[0] || data.model_name
          repos_map.set(repo_url, {
            name: label,
            model_key: data.model_key,
            repo: repo_url,
          })
        }
      } catch (error) {
        console.warn(
          `  Warning: ${yml_file} - ${error instanceof Error ? error.message : error}`,
        )
      }
    }
  }

  return Array.from(repos_map.values())
}

const is_cache_fresh = async (cache_file: string): Promise<boolean> => {
  if (force_refresh) return false
  try {
    return Date.now() - (await stat(cache_file)).mtimeMs < CACHE_TTL_MS
  } catch {
    return false
  }
}

const get_cached = async (cache_file: string): Promise<RepoData | null> => {
  if (!(await is_cache_fresh(cache_file))) return null
  try {
    const cached = JSON.parse(await readFile(cache_file, `utf-8`))
    return cached.name && cached.repo && cached.stars !== undefined ? cached : null
  } catch {
    return null
  }
}

const fetch_github = async (url: string, headers: Record<string, string>) => {
  const res = await fetch(url, { headers })
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  return res
}

const get_count_from_pagination = (link: string | null) => {
  const match = /page=(?<page>\d+)>; rel="last"/.exec(link ?? ``)
  return match?.groups ? Math.trunc(Number(match.groups.page)) : null
}

// total count from the Link header's last page, falling back to the response body length
const count_from_response = async (res: Response): Promise<number> =>
  get_count_from_pagination(res.headers.get(`Link`)) ?? (await res.json()).length ?? 0

const get_github_stats = async (
  model_info: ModelInfo,
  token: string | null,
  cache_file: string,
): Promise<RepoData | null> => {
  const { name, model_key, repo } = model_info
  const cached = await get_cached(cache_file)
  if (cached) {
    // Update model_key in case it changed (cache only stores GitHub stats)
    cached.model_key = model_key
    cached.name = name
    console.info(`✓ ${repo} (cached)`)
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
    const commits_last_year = await count_from_response(commits_res)

    // Fetch contributor count
    const contrib_res = await fetch_github(
      `https://api.github.com/repos/${repo}/contributors?per_page=1&anon=true`,
      headers,
    )
    const contributors = await count_from_response(contrib_res)
    const stars = repo_data.stargazers_count ?? 0
    const forks = repo_data.forks_count ?? 0
    const data: RepoData = {
      name,
      model_key,
      repo,
      stars,
      forks,
      commits_last_year,
      contributors,
    }

    // Cache the result
    await mkdir(cache_dir, { recursive: true })
    await writeFile(cache_file, JSON.stringify(data))
    console.info(
      `✓ ${repo}: ${data.stars}★ ${data.forks}f ${data.commits_last_year}c ${data.contributors}contrib`,
    )
    return data
  } catch (error) {
    console.error(`✗ ${repo}:`, error instanceof Error ? error.message : error)
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

  for (const info of repos) {
    const { repo } = info
    const existing = existing_map.get(repo)
    const cache_file = join(cache_dir, `${repo.replace(`/`, `_`)}.json`)
    const merged_existing = existing ? { ...existing, ...info } : null

    if (merged_existing && (await is_cache_fresh(cache_file))) {
      results.push(merged_existing)
      continue
    }

    const stats = await get_github_stats(info, token, cache_file)
    const result = stats ?? merged_existing
    if (result) results.push(result)
    await new Promise((resolve) => setTimeout(resolve, 300))
  }

  await writeFile(output_file, `${JSON.stringify(results, null, 2)}\n`)
  console.info(`Saved ${results.length}/${repos.length} repos`)
}

main().catch((error) => {
  console.error(`Fatal:`, error)
  process.exitCode = 1
})
