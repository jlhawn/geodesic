import { Figure, slider, choice, legend, caption, readout, text, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { loadColumn, K } from '../radiativeColumn.module.js';
import { SOLAR_CONSTANT } from '../../js/physics/radiation.module.js';
import { createSplit } from '../sunlightSplit.module.js';

const ZOOM = 10, TOP = 30, RATE = 30, PAD = 12, BARS = 184, WIDE = 560, GAP = 2, HALO = 'rgba(20, 20, 22, 0.85)';
const SKY = '#dfe6ee', SHINE = '#9aa3ad', SOIL = '#80664a', EARTH = '#5a4a36';
const PRESSURES = [1000, 300, 100, 30, 10, 3, 1, 0.3], HEIGHTS = [10, 20, 30, 40, 50];
const GROUNDS = [['sea', 'sea'], ['grassland', 'grassland'], ['forest', 'forest'], ['bare soil', 'soil'], ['fresh snow', 'snow'], ['sea ice', 'ice']];
const PARTS = ['byAir', 'byGround', 'ozone', 'gases', 'aerosol', 'ground'];

const mix = (hex, to, t) => `rgb(${[1, 3, 5].map((n) => Math.round(parseInt(hex.slice(n, n + 2), 16) * (1 - t) + to * t)).join(', ')})`;
const number = (v, digits = 1) => v.toLocaleString('en-US', { minimumFractionDigits: digits, maximumFractionDigits: digits });
const share = (v, of) => `${(100 * v / of).toFixed(1)}%`;
const hPa = (p) => { const v = p / 100; return v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v.toFixed(2); };

function segments(ctx, x, y, height, per, parts) {
  for (const [value, color] of parts) {
    const width = Math.max(0, value) * per;
    if (value > 0) {
      const inset = width > 2 * GAP + 2 ? GAP / 2 : 0, w = Math.max(1, width - 2 * inset);
      ctx.fillStyle = color;
      ctx.beginPath(); ctx.roundRect(x + inset, y, w, height, Math.min(3, w / 2)); ctx.fill();
    }
    x += width;
  }
}

function place(ctx, items, lo, hi, size) {
  ctx.font = `400 ${size}px system-ui, sans-serif`;
  for (const item of items) item.width = Math.max(ctx.measureText(item.name).width, item.value ? ctx.measureText(item.value).width : 0);
  let edge = lo;
  for (const item of items) { item.left = Math.max(edge, item.anchor - item.width / 2); edge = item.left + item.width + 10; }
  edge = hi;
  for (const item of [...items].reverse()) { item.left = Math.min(item.left, edge - item.width); edge = item.left - 10; }
  for (const item of items) item.center = Math.max(lo, item.left) + item.width / 2;
}

function labels(ctx, items, lo, hi, size, from, to, first) {
  place(ctx, items, lo, hi, size);
  const line = size + 2.5;
  for (const item of items) {
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(item.anchor, from); ctx.lineTo(item.center, to); ctx.stroke();
    const rows = item.value ? [[item.name, MUTED], [item.value, INK]] : [[item.name, MUTED]], start = to < from ? first + (2 - rows.length) * line : first;
    rows.forEach(([words, color], n) => text(ctx, words, item.center, start + n * line, { align: 'center', color, size }));
  }
}

export function mountBeam(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), GAS = mix(A, 255, 0.45), HAZE = mix(A, 0, 0.35);
  const state = { elevation: 90, ground: 'sea', ozone: true };
  let split = null, goal = null, shown = null, probe = null, plot = null;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 500 : 480), step, draw });
  caption(root, 'One column of dry air, the standard midlatitude summer with its water vapor taken out, lit through the model’s own radiation code. The top bar splits the sunlight arriving at the top of the air by where it ends up, on a track as long as the 1,362 W/m² of a sun straight overhead, and the air’s share is drawn again below it, ten times larger. Underneath is how fast the sunlight warms each of the model’s layers, in K per day, while the sun stands at that height; the topmost layer, above 0.3 hPa, is left off. Point at or tap the profile to read another layer.');
  legend(root, [['bar', 'absorbed in the air', A], ['bar', 'reflected by the air', SKY], ['bar', 'reflected by the ground', SHINE], ['bar', 'absorbed by the ground', SOIL], ['line', 'heating by sunlight', A], ['dash', 'tropopause', MUTED]]);
  slider(controls, { label: 'Sun’s height', min: 5, max: 90, step: 1, value: state.elevation, format: (v) => `${v}°, μ = ${Math.sin(v * Math.PI / 180).toFixed(2)}`, onInput: (v) => { state.elevation = v; update(); } });
  choice(controls, { label: 'Ground', options: GROUNDS, value: state.ground, onChange: (v) => { state.ground = v; update(); } });
  choice(controls, { label: 'Ozone', options: [['on', 'on'], ['off', 'off']], value: 'on', onChange: (v) => { state.ozone = v === 'on'; update(); } });
  const out = readout(controls);

  function update() {
    if (!split) return;
    const s = split(state);
    goal = { ...s, parts: Float64Array.from(PARTS, (name) => s[name]), heating: Float64Array.from(s.heating) };
    if (!shown) shown = { parts: Float64Array.from(goal.parts), heating: Float64Array.from(goal.heating) };
    const g = goal;
    out.set([
      ['arriving', `${number(g.arriving, 0)} W/m²`],
      ['reflected to space', `${number(g.reflected)} W/m², ${share(g.reflected, g.arriving)}`],
      ['absorbed in the air', `${number(g.air)} W/m², ${share(g.air, g.arriving)}`],
      ['absorbed by the ground', `${number(g.ground)} W/m², ${share(g.ground, g.arriving)}`],
      ['the ground’s albedo', g.albedo[0] === g.albedo[1] ? g.albedo[0].toFixed(2) : `${g.albedo[0].toFixed(3)} for the direct beam, ${g.albedo[1].toFixed(2)} for diffuse light`],
      ['peak heating', `${g.peak.rate.toFixed(g.peak.rate < 1 ? 2 : 1)} K/day at ${hPa(g.peak.pressure)} hPa, ${(g.peak.height / 1000).toFixed(0)} km up`],
      ['in the lowest 10 km', `${g.lower.low.toFixed(2)} to ${g.lower.high.toFixed(2)} K/day`],
    ]);
    fig.play(true);
    fig.render();
  }

  function step(dt) {
    if (!goal) return false;
    const f = 1 - Math.exp(-dt * 12);
    let moved = 0;
    for (const key of ['parts', 'heating']) for (let n = 0; n < goal[key].length; n++) {
      moved = Math.max(moved, Math.abs(goal[key][n] - shown[key][n]));
      shown[key][n] += (goal[key][n] - shown[key][n]) * f;
    }
    if (moved < 1e-3) { shown.parts.set(goal.parts); shown.heating.set(goal.heating); fig.play(false); }
    return true;
  }

  function bars(ctx, w, size) {
    const x0 = PAD, span = w - 2 * PAD, per = span / SOLAR_CONSTANT, zoom = per * ZOOM, yBar = 66, yZoom = 118, bh = 18;
    const [byAir, byGround, ozone, gases, aerosol, ground] = shown.parts, air = ozone + gases + aerosol, mine = air + byAir;
    text(ctx, `sunlight arriving at the top: ${number(goal.arriving, 0)} W/m²`, x0, 12, { color: INK, size: size + 0.5, weight: 500 });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.05)'; ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.roundRect(x0 + 0.5, yBar + 0.5, span - 1, bh - 1, 3); ctx.fill(); ctx.stroke();
    segments(ctx, x0, yBar, bh, per, [[air, A], [byAir, SKY], [byGround, SHINE], [ground, SOIL]]);
    const xAir = x0 + mine * per, xShine = xAir + byGround * per, xEnd = xShine + ground * per;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0, yBar - 3); ctx.lineTo(x0, yBar - 6); ctx.lineTo(xAir, yBar - 6); ctx.lineTo(xAir, yBar - 3); ctx.stroke();
    labels(ctx, [
      { anchor: (x0 + xAir) / 2, name: 'the air' },
      { anchor: (xAir + xShine) / 2, name: 'reflected by the ground', value: number(goal.byGround) },
      { anchor: (xShine + xEnd) / 2, name: 'absorbed by the ground', value: number(goal.ground) },
    ], x0, x0 + span, size, yBar - 7, yBar - 13, yBar - 34);
    const reach = Math.min(span, mine * zoom), squeeze = mine > 0 ? reach / (mine * zoom) : 1, yGap = yBar + bh, mid = (yGap + yZoom) / 2;
    ctx.fillStyle = 'rgba(255, 255, 255, 0.04)'; ctx.beginPath(); ctx.moveTo(x0, yGap); ctx.lineTo(xAir, yGap); ctx.lineTo(x0 + reach, yZoom); ctx.lineTo(x0, yZoom); ctx.closePath(); ctx.fill();
    ctx.strokeStyle = LINE; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(xAir, yGap + 2); ctx.lineTo(x0 + reach, yZoom - 2); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, 'ten times larger', x0 + span, mid, { align: 'right', color: MUTED, size });
    segments(ctx, x0, yZoom, bh, zoom * squeeze, [[ozone, A], [gases, GAS], [aerosol, HAZE], [byAir, SKY]]);
    let at = x0;
    const items = [['ozone', ozone, goal.ozone], ['O₂ + CO₂', gases, goal.gases], ['aerosol', aerosol, goal.aerosol], ['reflected by the air', byAir, goal.byAir]].map(([name, value, exact]) => {
      const width = Math.max(0, value) * zoom * squeeze, item = { anchor: at + width / 2, name, value: number(exact) };
      at += width;
      return item;
    });
    labels(ctx, items, x0, x0 + span, size, yZoom + bh + 1, yZoom + bh + 8, yZoom + bh + 15);
  }

  function profile(ctx, w, h, size) {
    const narrow = w < WIDE, left = narrow ? 56 : 62, right = w - (narrow ? 40 : 46), top = BARS + 30, bottom = h - 34;
    const ip = goal.interfacePressure, ps = ip[K], span = Math.log(ps) - Math.log(TOP);
    const y = (p) => top + (Math.log(p) - Math.log(TOP)) / span * (bottom - top), x = (r) => left + r / RATE * (right - left);
    const at = (z) => { const { height: zs, pressure: p } = goal; if (z <= zs[K - 1]) return ps * Math.exp(-z / zs[K - 1] * Math.log(ps / p[K - 1])); for (let k = K - 1; k > 0; k--) if (zs[k - 1] >= z) return Math.exp(Math.log(p[k]) + (z - zs[k]) / (zs[k - 1] - zs[k]) * (Math.log(p[k - 1]) - Math.log(p[k]))); return p[0]; };
    plot = { y, top, bottom, left, right };
    text(ctx, 'heating by sunlight, K per day', left, BARS + 10, { color: INK, size: size + 0.5, weight: 500 });
    ctx.lineWidth = 1;
    for (let r = 0; r <= RATE; r += 10) { ctx.strokeStyle = r ? GRID : LINE; ctx.beginPath(); ctx.moveTo(x(r), top); ctx.lineTo(x(r), bottom); ctx.stroke(); text(ctx, `${r}`, x(r), bottom + 16, { align: 'center', color: MUTED, size }); }
    for (const p of PRESSURES) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(p * 100)); ctx.lineTo(right, y(p * 100)); ctx.stroke(); text(ctx, `${p} hPa`, left - 6, y(p * 100), { align: 'right', color: MUTED, size }); }
    for (const z of HEIGHTS) { const yy = y(at(z * 1000)); ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(right, yy); ctx.lineTo(right + 4, yy); ctx.stroke(); text(ctx, `${z} km`, right + 7, yy, { color: MUTED, size }); }
    ctx.fillStyle = EARTH; ctx.fillRect(left, bottom, right - left, 5);
    if (goal.tropopause) {
      const yy = y(goal.tropopause.pressure);
      ctx.strokeStyle = MUTED; ctx.setLineDash([3, 4]); ctx.beginPath(); ctx.moveTo(left, yy); ctx.lineTo(right, yy); ctx.stroke(); ctx.setLineDash([]);
      text(ctx, 'tropopause', right - 4, yy - 8, { align: 'right', color: MUTED, size, halo: HALO });
    }
    const r = shown.heating;
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    ctx.beginPath(); ctx.moveTo(x(0), y(ip[1]));
    for (let k = 1; k < K; k++) { ctx.lineTo(x(r[k]), y(ip[k])); ctx.lineTo(x(r[k]), y(ip[k + 1])); }
    ctx.lineTo(x(0), y(ps)); ctx.closePath();
    ctx.globalAlpha = 0.18; ctx.fillStyle = A; ctx.fill(); ctx.globalAlpha = 1;
    ctx.beginPath(); ctx.moveTo(x(r[1]), y(ip[1]));
    for (let k = 1; k < K; k++) { ctx.lineTo(x(r[k]), y(ip[k])); ctx.lineTo(x(r[k]), y(ip[k + 1])); }
    ctx.strokeStyle = A; ctx.lineWidth = 2; ctx.lineJoin = 'round'; ctx.stroke();
    ctx.restore();
    const low = goal.lower, yLow = y(at(4000));
    text(ctx, `${low.low.toFixed(2)} to ${low.high.toFixed(2)} K per day in the lowest 10 km`, x(0) + 8, yLow, { color: MUTED, size, halo: HALO });
    const mark = (k, strong) => {
      const px = x(r[k]), py = y(goal.pressure[k]), rate = `${goal.heating[k].toFixed(goal.heating[k] < 1 ? 2 : 1)} K/day at ${hPa(goal.pressure[k])} hPa`;
      ctx.font = `400 ${size}px system-ui, sans-serif`;
      const options = [`${rate}, ${(goal.height[k] / 1000).toFixed(0)} km up`, rate].flatMap((label) => [[label, px + 10, 'left'], [label, px - 10, 'right']]);
      const fits = ([label, from, align]) => { const width = ctx.measureText(label).width; return align === 'left' ? from + width <= right - 2 : from - width >= left + 2; };
      const [label, from, align] = options.find(fits) ?? [rate, clamp(px - ctx.measureText(rate).width / 2, left + 2, right - 2 - ctx.measureText(rate).width), 'left'];
      if (!strong) { ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]); }
      ctx.fillStyle = strong ? A : INK; ctx.strokeStyle = HALO; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2); ctx.stroke(); ctx.fill();
      text(ctx, label, from, options.some(fits) ? py : py + 14, { align, color: INK, size, halo: HALO });
    };
    if (probe !== null && probe !== goal.peak.k) mark(probe, false);
    mark(goal.peak.k, true);
  }

  function draw(ctx, w, h) {
    const size = w < WIDE ? 10.5 : 11;
    if (!shown) { text(ctx, 'loading the standard atmosphere', w / 2, h / 2, { align: 'center', color: MUTED, size }); return; }
    bars(ctx, w, size);
    profile(ctx, w, h, size);
  }

  const pick = (event) => {
    if (!plot || !goal) return;
    const r = fig.canvas.getBoundingClientRect(), py = event.clientY - r.top, px = event.clientX - r.left;
    if (py < plot.top - 8 || py > plot.bottom + 8 || px < plot.left - 40) { if (probe !== null) { probe = null; fig.render(); } return; }
    let best = 1;
    for (let k = 1; k < K; k++) if (Math.abs(plot.y(goal.pressure[k]) - py) < Math.abs(plot.y(goal.pressure[best]) - py)) best = k;
    if (best !== probe) { probe = best; fig.render(); }
  };
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (event) => { if (event.pointerType === 'mouse') pick(event); });
  fig.canvas.addEventListener('pointerleave', (event) => { if (event.pointerType === 'mouse' && probe !== null) { probe = null; fig.render(); } });

  loadColumn().then((profile) => { split = createSplit(profile); update(); });
}
