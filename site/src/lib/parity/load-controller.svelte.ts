import { get_error_message, type LoadStatus } from '$lib/asset-loader'

/** Shared race cancellation and status handling for asynchronous parity plots. */
export class ParityLoadController {
  status = $state<LoadStatus>(`idle`)
  error_message = $state(``)
  private load_id = 0

  /** Mark the current plot data as ready without starting a new request. */
  set_ready = (): void => {
    this.load_id++
    this.status = `ready`
    this.error_message = ``
  }

  /** Cancel pending requests and expose a terminal load error. */
  set_error = (message: string): void => {
    this.load_id++
    this.status = `error`
    this.error_message = message
  }

  /** Run one load, ignoring completion or errors after a newer load starts. */
  run = async (load: (is_current: () => boolean) => Promise<void>): Promise<void> => {
    const current_load_id = ++this.load_id
    const is_current = (): boolean => current_load_id === this.load_id
    this.status = `loading`
    this.error_message = ``
    try {
      await load(is_current)
      if (is_current()) this.status = `ready`
    } catch (error) {
      if (!is_current()) return
      this.status = `error`
      this.error_message = get_error_message(error)
    }
  }
}
