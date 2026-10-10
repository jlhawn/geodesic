import { mountCrossSection } from './crossSection.module.js';
import { jw06Section, heldSuarezEquilibrium, zonalMeanSection } from '../layerCases.module.js';
import { loadFrames } from '../frames.module.js';

const SHADE = ['ramp', 'temperature, 200 to 310 K', 'rgb(52, 55, 62)', 'warm'];
const ISENTROPES = ['faint', 'potential temperature every 10 K', 'rgba(255,255,255,0.6)'];
const WIND = ['line', 'eastward wind every 5 m/s, westward dashed', 'rgb(120, 220, 255)'];

export function mountJwSection(root) {
  const data = jw06Section();
  let top = -Infinity, at = 0;
  data.wind.forEach((column, i) => column.forEach((u, j) => { if (u > top) { top = u; at = j; } }));
  mountCrossSection(root, {
    modes: [{ label: 'the jets', data: () => data, readout: () => [['fastest wind', `${top.toFixed(0)} m/s near ${data.pressures[at].toFixed(0)} hPa`]] }],
    legendFor: () => [SHADE, ISENTROPES, WIND],
  });
}

export function mountHsSection(root) {
  const equilibrium = heldSuarezEquilibrium();
  let settled = null;
  const view = mountCrossSection(root, {
    modes: [
      { label: 'what the sunlight asks for', data: () => equilibrium, readout: () => [['no wind', 'this temperature pattern has no motion of its own']] },
      { label: 'what the atmosphere settles on', data: () => settled, readout: () => settled ? summary(settled) : [['loading', '…']] },
    ],
    initial: 0,
    legendFor: (mode) => (mode ? [SHADE, ISENTROPES, WIND] : [SHADE, ISENTROPES]),
  });
  loadFrames(new URL('../data/heldsuarez_N16.bin', import.meta.url).href).then((frames) => { settled = zonalMeanSection(frames.header.zonalMean); view.refresh(); });
}

function summary(d) {
  let best = { u: -Infinity }, trades = Infinity;
  d.wind.forEach((column, i) => column.forEach((u, j) => { if (u > best.u) best = { u, lat: d.lats[i], p: d.pressures[j] }; }));
  d.wind.forEach((column, i) => { if (Math.abs(d.lats[i]) < 20) trades = Math.min(trades, column[column.length - 1]); });
  return [['strongest jet', `${best.u.toFixed(0)} m/s at ${Math.abs(best.lat).toFixed(0)}°${best.lat > 0 ? 'N' : 'S'}, ${best.p.toFixed(0)} hPa`], ['surface wind in the tropics', `${trades.toFixed(1)} m/s`]];
}
