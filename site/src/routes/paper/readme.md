# Edit the paper interactively

Install [`npm`](https://npmjs.com), then from the `site/` directory run

```sh
npm install
npm run dev
```

Both steps should only take a few seconds after which the terminal should show something like

```text
  VITE v4.0.4  ready in 531 ms

  ➜  Local:   http://localhost:3000/
  ➜  Network: use --host to expose
  ➜  press h to show help
```

Open the development server running on <http://localhost:3000/paper> to view the paper and see hot-reloaded changes whenever you edit and save the file.

## Auto-number headings

Currently done in CSS (see `headings-number.css`) but could be done with JS using

```js
<script>
  import { onMount } from 'svelte'

  onMount(() => {
    const counter = { h2: 0, h3: 0, h4: 0 }
    // TODO define last_level to reset h3 counter on next h2 heading

    for (const heading of document.querySelectorAll(`h2, h3, h4`)) {
      const level = heading.tagName.toLowerCase()
      counter[level]++
      const { h2, h3, h4 } = counter

      if (level === `h2`) {
        heading.innerHTML = `${h2} ${heading.innerHTML}`
      } else if (level === `h3`) {
        heading.innerHTML = `${h2}.${h3} ${heading.innerHTML}`
      } else if (level === `h4`) {
        heading.innerHTML = `${h2}.${h3}.${h4} ${heading.innerHTML}`
      } else throw `Unexpected heading tag ${level}`
    }
  })
</script>
```

## Convert Markdown to LaTeX

Install [`pandoc`](https://pandoc.org/installing.html) and run

```sh
cd site/src/routes/paper
pandoc -s +page.md -o paper.tex --wrap=preserve
```

## Markdown in LaTeX

A better option than converting back and forth between TeX and Markdown might be to use [this Overleaf template](https://overleaf.com/latex/examples/using-markdown-in-latex-documents/whdrnpcpnwrm) which enables importing Markdown files into TeX documents.
