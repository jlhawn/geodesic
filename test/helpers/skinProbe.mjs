export const SEEN = 5, STEP = 8;

const views = ({ seen, seenSet, stepped, steppedSet, flux }) => ({ seen: new Float64Array(seen), seenSet: new Uint8Array(seenSet), stepped: new Float64Array(stepped), steppedSet: new Uint8Array(steppedSet), flux: new Float64Array(flux) });

export function skinRecords(C) {
  const data = { seen: new SharedArrayBuffer(8 * SEEN * C), seenSet: new SharedArrayBuffer(C), stepped: new SharedArrayBuffer(8 * STEP * C), steppedSet: new SharedArrayBuffer(C), flux: new SharedArrayBuffer(8 * C) };
  return { data, ...views(data) };
}

export default function probe(model, data) {
  const { seen, seenSet, stepped, steppedSet, flux } = views(data);
  const { exchange, skin, state, radiation, phases } = model, original = exchange.cell, step = skin.step, physics = phases.physics;
  exchange.cell = (i, pi, theta, q, qc, surfaceT, lowest, concentration, ...rest) => {
    if (model.geography.land[i] || concentration > 0 || state[6][i] > 0) { seenSet[i] = 0; return original(i, pi, theta, q, qc, surfaceT, lowest, concentration, ...rest); }
    const slab = original(i, pi, theta, q, qc, state[3][i], lowest, concentration, ...rest).ri;
    const out = original(i, pi, theta, q, qc, surfaceT, lowest, concentration, ...rest);
    seen.set([slab, out.ri, out.z, out.momentum, out.heat], SEEN * i);
    seenSet[i] = 1;
    return out;
  };
  skin.step = (i, ...args) => { stepped[STEP * i] = skin.warm[i]; stepped.set(args, STEP * i + 1); steppedSet[i] = 1; step(i, ...args); };
  phases.physics = (iFrom, iTo, ...rest) => { physics(iFrom, iTo, ...rest); flux.set(radiation.surfaceFlux.subarray(iFrom, iTo), iFrom); };
}
