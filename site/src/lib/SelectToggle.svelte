<script lang="ts">
  import { Icon } from 'svelte-widgets'
  import { Spinner } from 'matterviz/feedback'
  import { tooltip as tip } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  interface OptionInfo {
    value: string
    label: string
    tooltip?: string
    link?: string
    loading?: boolean // show a small spinner in the button (e.g. while its tab loads)
  }
  let {
    selected = $bindable(``),
    options = [],
    tooltip_placement = `bottom`,
    ...rest
  }: {
    selected: string
    options: OptionInfo[]
    tooltip_placement?: `top` | `bottom` | `left` | `right`
  } & HTMLAttributes<HTMLDivElement> = $props()
  const target = { target: `_blank`, rel: `noopener noreferrer` }
</script>

<div role="group" {...rest} class={[`selection-toggle`, rest.class]}>
  {#each options as { value, label, tooltip, link, loading } (value)}
    <!-- the info link is a sibling of the button, not a child: interactive content
    inside a <button> is invalid HTML. .option carries the pill styling so the icon
    still renders inside the border. -->
    <span class="option" class:active={selected === value}>
      <button
        aria-pressed={selected === value}
        onclick={() => (selected = value)}
        {@attach tip({
          allow_html: true,
          content: tooltip,
          placement: tooltip_placement,
        })}
      >
        {@html label}
        {#if loading}
          <Spinner
            style="--spinner-size: 0.9em; --spinner-border-width: 2px; --spinner-margin: 0"
          />
        {/if}
      </button>
      {#if link}
        <a href={link} {...target}>
          <Icon icon="Info" style="transform: translateY(-1px)" />
        </a>
      {/if}
    </span>
  {/each}
</div>

<style>
  .selection-toggle {
    display: flex;
    flex-wrap: wrap;
    place-content: center;
    gap: 8pt;
  }
  /* the pill: border and state colors live here so the info link renders inside them,
  but padding stays on the button so clicking the pill's edges still toggles */
  .selection-toggle .option {
    display: inline-flex;
    align-items: stretch;
    border: 0.5px solid var(--border);
    border-radius: 3pt;
    color: var(--btn-text);
  }
  .selection-toggle .option:hover {
    background: var(--nav-bg);
  }
  .selection-toggle .option.active {
    border-color: var(--link-color);
    background: color-mix(in oklab, var(--link-color) 8%, var(--nav-bg));
    color: var(--link-color);
  }
  .selection-toggle button,
  .selection-toggle button:hover {
    display: inline-flex;
    align-items: center;
    gap: 0.5ex;
    padding: 4px 8px;
    border: none;
    border-radius: 0;
    background: transparent;
    color: inherit;
  }
  /* 0.83em ≈ the UA font size <button> labels render at, so the info icon (which now
  sits outside the button) still matches them instead of the larger page font. The
  button's right padding shrinks to the gap the icon used to have inside it. */
  .selection-toggle .option:has(a) button {
    padding-right: 0.5ex;
  }
  .selection-toggle .option a {
    display: inline-flex;
    align-items: center;
    padding-right: 8px;
    font-size: 0.83em;
  }
</style>
