<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'

  export let show_tips: boolean = false
  export let tips_title = `Usage Tips`
  export let dialog: HTMLDialogElement
  export let btn: HTMLButtonElement

  function close_if_outside_click(event: MouseEvent) {
    const is_outside = dialog && !dialog.contains(event.target as Node)
    if (show_tips && is_outside && !btn.contains(event.target as Node)) {
      show_tips = false
    }
  }

  function toggle(event: KeyboardEvent) {
    if (event.key == `Escape`) show_tips = false
    if (event.key == `j` && event.metaKey) show_tips = !show_tips
  }
</script>

<svelte:window on:click={close_if_outside_click} on:keydown={toggle} />

<footer>
  <nav>
    <a href="{repository}/issues">Issues</a>
    <a href="mailto:janosh@lbl.gov?subject=Matbench Discovery">Contact</a>
    <a href="/changelog">Changelog</a>
    <button
      on:click={() => (show_tips = true)}
      bind:this={btn}
      title={tips_title}
      style="padding: 0; transform: scale(1.2);"
    >
      <Icon icon="mdi:lightbulb-on-outline" inline />
    </button>
  </nav>
  <img
    src="/favicon.svg"
    alt="Logo"
    width="30px"
    style="vertical-align:middle;"
  />&emsp;Matbench Discovery
</footer>

<dialog bind:this={dialog} open={show_tips}>
  <h3>{tips_title}</h3>
  <p title="For keyboard-only site navigation">
    <kbd>cmd</kbd> + <kbd>k</kbd> to bring up a nav palette.
  </p>
  <p title="For keyboard-only site navigation">
    <kbd>cmd</kbd> + <kbd>j</kbd> to bring up these site options.
  </p>
  <p title="single click to hide a trace">
    <kbd>single click</kbd> a figure legend handle to hide its corresponding trace in a plot
  </p>
  <p title="double click to hide all other traces">
    <kbd>double click</kbd> the figure legend handle of any trace in a plot to hide all others.
    double click again to show all.
  </p>
</dialog>

<style>
  footer {
    padding: 3vh 3vw;
    background: #00061a;
    text-align: center;
  }
  footer nav {
    display: flex;
    gap: 1em;
    justify-content: center;
    margin: 2em 0;
  }
  footer > nav > button {
    background: none;
    color: var(--blue);
  }
  dialog {
    visibility: hidden;
    opacity: 0;
    position: fixed;
    top: 40%;
    background: var(--sms-options-bg);
    color: var(--text-color);
    display: grid;
    border: none;
    border-radius: 3pt;
    padding: 1ex 2ex;
    transition: 0.2s;
  }
  dialog[open] {
    visibility: visible;
    opacity: 1;
  }
  dialog > :is(p, h3) {
    margin: 0;
  }
</style>
