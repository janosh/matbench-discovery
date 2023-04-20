<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'

  let show_tips: boolean
  const tips_title = `Usage Tips`
  let dialog: HTMLDialogElement
  let btn: HTMLButtonElement

  function close_if_outside_click(event: MouseEvent) {
    if (
      dialog &&
      show_tips &&
      !dialog.contains(event.target as Node) &&
      !btn.contains(event.target as Node)
    ) {
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
  Questions/feedback?
  <a href="{repository}/issues"><Icon icon="octicon:mark-github" inline /></a>
  <a href="mailto:janosh@lbl.gov?subject=Matbench Discovery">
    <Icon icon="mdi:email" inline />
  </a>
  <!-- open modal on clicking tips icon -->
  <button on:click={() => (show_tips = true)} bind:this={btn} title={tips_title}>
    <Icon icon="mdi:lightbulb-on-outline" inline />
  </button>
</footer>

<dialog bind:this={dialog} open={show_tips}>
  <h3>{tips_title}</h3>
  <p title="For keyboard-only site navigation">
    Use <kbd>cmd+k</kbd> to bring up a nav palette.
  </p>
  <p title="For keyboard-only site navigation">
    Use <kbd>cmd+j</kbd> to bring up these site options.
  </p>
</dialog>

<style>
  footer {
    padding: 3vh 3vw;
    background: #00061a;
    text-align: center;
  }
  footer > a {
    margin: 0 0 0 4pt;
  }
  footer > button {
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
  dialog > * {
    margin: 0;
  }
  p kbd {
    background-color: rgba(255, 255, 255, 0.1);
    padding: 1pt 3pt;
    font-size: larger;
    border-radius: 2pt;
  }
</style>
