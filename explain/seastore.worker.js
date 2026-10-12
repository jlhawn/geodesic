import { createSeaStore } from './seaStore.module.js';

const BUDGET_MS = 30;
let store = null, run = 0, timer = 0;

const transferables = (batch) => Object.values(batch).filter((v) => ArrayBuffer.isView(v)).map((v) => v.buffer);

function work() {
  timer = 0;
  if (!store) return;
  const t0 = performance.now();
  while (!store.done && performance.now() - t0 < BUDGET_MS) store.run(1);
  const batch = store.batch();
  self.postMessage({ type: 'samples', run, done: store.done, ...batch }, transferables(batch));
  if (!store.done) timer = setTimeout(work, 0);
}

self.onmessage = ({ data }) => {
  if (data.type !== 'start') return;
  clearTimeout(timer);
  run = data.run;
  try {
    store = createSeaStore(data.settings);
    self.postMessage({ type: 'layout', run, ...store.layout });
    work();
  } catch (error) {
    store = null;
    self.postMessage({ type: 'error', run, message: String(error && error.message ? error.message : error) });
  }
};
