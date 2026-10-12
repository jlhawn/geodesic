import { CASES } from './oceanCases.module.js';

const BUDGET_MS = 40;
let sim = null, run = 0, owed = 0, endless = false;

function transferables(value, list = []) {
  if (ArrayBuffer.isView(value)) list.push(value.buffer);
  else if (value && typeof value === 'object') for (const part of Object.values(value)) transferables(part, list);
  return list;
}

function post(type, payload = {}) {
  self.postMessage({ type, run, ...payload }, transferables(payload));
}

function start({ kind, options = {}, run: id = 0 }) {
  run = id; owed = 0; endless = !!options.endless; sim = null;
  if (!CASES[kind]) throw new Error(`no ocean case called ${kind}; the cases are ${Object.keys(CASES).join(', ')}`);
  sim = CASES[kind](options);
  post('layout', sim.layout());
  post('frame', { ...sim.frame(), steps: 0, cost: 0 });
}

function advance(seconds) {
  owed += seconds;
  const t0 = performance.now();
  let steps = 0;
  for (;;) {
    const dt = sim.nextDt();
    if (owed < dt || (!endless && sim.time >= sim.end) || performance.now() - t0 >= BUDGET_MS) break;
    sim.step(); owed -= dt; steps++;
  }
  const cost = steps ? (performance.now() - t0) / steps : 0;
  owed = Math.min(owed, sim.nextDt());
  post('frame', { ...sim.frame(), steps, cost });
}

self.onmessage = ({ data }) => {
  try {
    if (data.type === 'start') start(data);
    else if (data.type === 'advance' && sim) advance(data.seconds);
    else if (data.type === 'stop') sim = null;
  } catch (error) {
    sim = null;
    post('error', { message: String(error && error.message ? error.message : error) });
  }
};
