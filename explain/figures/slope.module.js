import { Figure, slider, choice, buttons, legend, caption, readout, text, arrow, rampRGB, anomalyColor, termColor, clamp, DARK_NEUTRAL, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { R, G, heightOf } from '../physics.module.js';

const LEVELS = [1000, 850, 700, 500, 300, 200], KM = 1000, WIDTH = 2000 * KM, TOP = 15000, START = 400 * KM, DX = 10 * KM, SAMPLES = 40, SIX = 6 * 3600, SPEED = SIX / 4, STEEP = 5, HALO = 'rgba(20, 20, 22, 0.85)';
const MAX_PUSH = R * 20 * Math.log(1000 / 200) / WIDTH, MAX_WIND = MAX_PUSH * SIX, CORIOLIS = 2 * 7.2921e-5 * Math.SQRT1_2, TURNED = CORIOLIS * SIX / 2;
const round2 = (n) => { const k = 10 ** Math.max(0, Math.floor(Math.log10(n)) - 1); return (Math.round(n / k) * k).toLocaleString('en-US'); };

export function mountSlope(root) {
  const controls = root.querySelector('.controls'), blue = termColor('b');
  let contrast = 10, level = 4, steep = true, t = 0, geometry = null;
  const offset = (m) => contrast * (0.5 - m / WIDTH);
  const height = (p, m) => heightOf(p * 100) + R * offset(m) / G * Math.log(1000 / p);
  const push = (p, m) => -G * (height(p, m + DX) - height(p, m - DX)) / (2 * DX);
  const parcel = (p) => { const a = push(p, START); return { m: START + 0.5 * a * t * t, u: a * t }; };

  const fig = new Figure(root, { height: 400, minHeight: 340, step, draw });
  caption(root, 'A slice of atmosphere 2,000 km across, warm on the left and cold on the right, with the same 1000 hPa at the ground on both sides. The dots are parcels of air, one on each surface of equal pressure, at rest until you let them go.');
  const items = () => [['force', 'the push along the chosen surface: g times its downhill slope', blue], ['arrow', 'the wind it builds once the air is let go', INK], ['faint', steep ? `surfaces of equal pressure, tilted ${STEEP} times more than they really are` : 'surfaces of equal pressure, at their true heights', LINE], ['ramp', 'warmer to colder air', 'warm', 'cool', `rgb(${DARK_NEUTRAL.join(', ')})`]];
  const key = legend(root, items());
  slider(controls, { label: 'The warm side is warmer by', min: 0, max: 20, step: 1, value: contrast, format: (v) => `${v} °C`, onInput: (v) => { contrast = v; update(); } });
  const levelSlider = slider(controls, { label: 'Pressure surface', min: 0, max: LEVELS.length - 1, step: 1, value: level, format: (v) => (LEVELS[v] === 1000 ? '1000 hPa, at the ground' : `${LEVELS[v]} hPa, about ${(heightOf(LEVELS[v] * 100) / KM).toFixed(1)} km up`), onInput: (v) => { level = v; update(); } });
  choice(controls, { label: 'Tilt', options: [[`drawn ${STEEP} times steeper`, 'steep'], ['true to the height axis', 'true']], value: 'steep', onChange: (v) => { steep = v === 'steep'; key.set(items()); fig.render(); }, span: true });
  buttons(controls, [['Let go for six hours', () => { t = 0; fig.play(true); fig.render(); }], ['Reset', () => { t = 0; fig.play(false); fig.render(); }]]);
  const out = readout(controls);

  function update() { report(); fig.render(); }

  function report() {
    const p = LEVELS[level], drop = height(p, 0) - height(p, WIDTH), a = push(p, WIDTH / 2), name = p === 1000 ? 'the ground, 1000 hPa' : `the ${p} hPa surface`;
    if (drop < 0.05) return out.set([[name, 'the same height on both sides'], ['slope', 'none'], ['push', 'none'], ['after six hours', 'still calm']]);
    const hour = a * 3600, speed = (u) => `${u.toFixed(u < 10 ? 1 : 0)} m/s`;
    out.set([[name, `${drop.toFixed(0)} m lower on the cold side`], ['slope', `${(drop / WIDTH * 100 * KM).toFixed(1)} m per 100 km, 1 in ${round2(WIDTH / drop)}`], ['push toward the cold side, g times the slope', `${a.toPrecision(2)} m/s², ${hour.toFixed(hour < 1 ? 2 : 1)} m/s an hour`], ['wind after six hours if nothing else acted', speed(a * SIX)], ['with the Coriolis force at 45° north', `${speed(2 * a / CORIOLIS * Math.sin(TURNED))}, turned ${Math.round(TURNED * 180 / Math.PI)}° to its right`]]);
  }

  function step(dt) {
    if (t >= SIX) { fig.play(false); return false; }
    t = Math.min(SIX, t + dt * SPEED);
  }

  function draw(ctx, w, h) {
    const left = 40, right = 12, top = 24, bottom = h - 32, span = w - left - right, slots = w < 500 ? 4 : 6;
    const x = (m) => left + m / WIDTH * span, y = (z) => bottom - z / TOP * (bottom - top);
    const shown = (p, m) => { const mid = heightOf(p * 100); return y(mid + (steep ? STEEP : 1) * (height(p, m) - mid)); };
    geometry = { left, span, shown };
    const tint = ctx.createLinearGradient(left, 0, left + span, 0), rgb = [0, 0, 0];
    for (const f of [0, 0.5, 1]) { rampRGB(0.5 + offset(f * WIDTH) / 40, rgb, DARK_NEUTRAL); tint.addColorStop(f, `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`); }
    ctx.fillStyle = tint; ctx.fillRect(left, top, span, bottom - top);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.setLineDash([2, 4]);
    for (let z = 2000; z < TOP; z += 2000) { ctx.beginPath(); ctx.moveTo(left, y(z)); ctx.lineTo(left + span, y(z)); ctx.stroke(); }
    ctx.setLineDash([]);
    for (let z = 0; z < TOP; z += 2000) text(ctx, z + 2000 >= TOP ? `${z / KM} km` : `${z / KM}`, left - 6, y(z), { align: 'right', color: MUTED, size: 10 });
    for (let d = 0; d <= 2000; d += 500) text(ctx, d === 2000 ? '2,000 km' : d.toLocaleString('en-US'), x(d * KM), bottom + 20, { align: d === 0 ? 'left' : d === 2000 ? 'right' : 'center', color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, span, 6);
    LEVELS.forEach((p, k) => {
      const chosen = k === level;
      ctx.strokeStyle = chosen ? INK : LINE; ctx.lineWidth = chosen ? 2 : 1.2;
      ctx.beginPath(); for (let s = 0; s <= SAMPLES; s++) { const m = s / SAMPLES * WIDTH; if (s) ctx.lineTo(x(m), shown(p, m)); else ctx.moveTo(x(m), shown(p, m)); } ctx.stroke();
      text(ctx, `${p} hPa`, left + span - 4, shown(p, WIDTH) - 9, { align: 'right', color: chosen ? INK : MUTED, size: 10, weight: chosen ? 600 : 400, halo: HALO });
    });
    const p = LEVELS[level], room = span / slots;
    for (let i = 0; i < slots; i++) {
      const m = (i + 0.5) / slots * WIDTH, len = push(p, m) / MAX_PUSH * 0.8 * room, cx = x(m), cy = Math.max(shown(p, m - len / 2 / span * WIDTH), shown(p, m + len / 2 / span * WIDTH)) + 9;
      if (len >= 1) arrow(ctx, cx - len / 2, cy, cx + len / 2, cy, { color: blue, width: 2, head: clamp(len * 0.6, 4, 8), dash: [5, 4], open: true });
    }
    const dots = LEVELS.map((q) => { const { m, u } = parcel(q); return { u, px: x(m), py: shown(q, m) }; });
    ctx.strokeStyle = 'rgba(255,255,255,0.25)'; ctx.lineWidth = 1;
    ctx.beginPath(); dots.forEach(({ px, py }, k) => { if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }); ctx.stroke();
    dots.forEach(({ u, px, py }, k) => {
      const chosen = k === level, len = u / MAX_WIND * 0.25 * span;
      arrow(ctx, px, py, px + len, py, { color: chosen ? INK : 'rgba(221,221,221,0.55)', width: chosen ? 2 : 1.5, head: chosen ? 8 : 6 });
      ctx.fillStyle = chosen ? '#fff' : 'rgba(255,255,255,0.6)'; ctx.beginPath(); ctx.arc(px, py, chosen ? 4.5 : 3, 0, Math.PI * 2); ctx.fill();
      if (chosen && t > 0) text(ctx, `${u.toFixed(u < 0.05 || u >= 10 ? 0 : 1)} m/s`, px + 6, py - 12, { color: INK, size: 10, halo: HALO });
    });
    if (contrast) { text(ctx, 'warm side', left + 2, 12, { color: anomalyColor(1), size: 11 }); text(ctx, 'cold side', left + span - 2, 12, { align: 'right', color: anomalyColor(-1), size: 11 }); }
    text(ctx, t > 0 ? `after ${(t / 3600).toFixed(1)} hours` : 'tap a line to choose it', left + span / 2, 12, { align: 'center', color: MUTED, size: 11 });
  }

  function pick({ x: px, y: py }, reach) {
    if (!geometry || px < geometry.left - 12 || px > geometry.left + geometry.span + 12) return null;
    const m = clamp((px - geometry.left) / geometry.span, 0, 1) * WIDTH;
    let best = null, gap = reach;
    LEVELS.forEach((p, k) => { const d = Math.abs(geometry.shown(p, m) - py); if (d < gap) { gap = d; best = k; } });
    return best;
  }
  const choose = (at, reach = Infinity) => { const k = pick(at, reach); if (k === null || k === level) return; level = k; levelSlider.value = k; update(); };
  fig.pointer({ hit: (at) => pick(at, 18) !== null, down: (at) => choose(at, 18), move: (at) => choose(at) });

  report();
}
