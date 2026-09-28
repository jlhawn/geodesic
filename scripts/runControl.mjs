import { spawn } from 'node:child_process';
import { existsSync } from 'node:fs';

/*
 * What a spin-up script needs to survive a preemptible machine.
 *
 * stopOnSignal turns the first SIGTERM or SIGINT into a request that the
 * stepping loop polls (`requested` holds the signal's name), so that the
 * script can finish its ocean step, save and exit 0. Later signals are
 * logged and otherwise ignored: a driver forwards the signal its process
 * group may already have delivered.
 *
 * syncAfterSave runs `command` (SYNC_CMD) through /bin/sh after each file
 * it is handed, with the file's path as $1, one at a time in the order
 * handed, each in its own process group so that a signal meant for the
 * spin-up does not cut an upload short. A file pruned before its turn is
 * skipped, since a newer one follows it in the queue. A failing command
 * is tried three times, 5 s apart, and then logged; `drain` resolves once
 * the queue is empty.
 */
export function stopOnSignal(log) {
  const control = { requested: null };
  for (const signal of ['SIGTERM', 'SIGINT']) {
    process.on(signal, () => {
      if (control.requested) { log(`${signal} while stopping; still saving`); return; }
      control.requested = signal;
      log(`${signal} at ${new Date().toISOString()}: saving after the ocean step in progress`);
    });
  }
  return control;
}

function runOnce(command, path) {
  return new Promise((resolve) => {
    const child = spawn('/bin/sh', ['-c', command, 'sync', path], { detached: true, stdio: ['ignore', 'ignore', 'pipe'] });
    let errors = '';
    child.stderr.on('data', (chunk) => { errors = (errors + chunk).slice(-400); });
    child.on('error', (error) => resolve({ code: -1, errors: error.message }));
    child.on('close', (code, signal) => resolve({ code: code ?? signal, errors }));
  });
}

export function syncAfterSave(command, log) {
  let queue = Promise.resolve();
  async function sync(path) {
    for (let attempt = 1; ; attempt++) {
      if (!existsSync(path)) return;
      const { code, errors } = await runOnce(command, path);
      if (code === 0) return;
      if (attempt === 3) { log(`SYNC_CMD failed on ${path} (exit ${code}) three times: ${errors.trim().replace(/\s+/g, ' ')}`); return; }
      await new Promise((resolve) => setTimeout(resolve, 5000));
    }
  }
  return {
    after(path) { if (command) queue = queue.then(() => sync(path)); },
    drain() { return queue; },
  };
}
