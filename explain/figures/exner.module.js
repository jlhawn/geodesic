import { Figure, slider, buttons, legend, caption, readout, text, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { R, CP, KAPPA, G, P0, T0, LAPSE, TROPOPAUSE, T_STRATOSPHERE } from '../physics.module.js';

const LEVELS = [1000, 850, 700, 500, 300, 200, 100], K = LEVELS.length - 1, RANGE = 20, ZMAX = 18000, EDGE = 0.5, GAIN = 1 / 3, SCALES = [250, 500, 1000], HALO = 'rgba(20, 20, 22, 0.85)';
const SEA = 101325, P11 = SEA * (T_STRATOSPHERE / T0) ** (G / (R * LAPSE));
const standard = (p) => p >= P11 ? T0 / LAPSE * (1 - (p / SEA) ** (R * LAPSE / G)) : TROPOPAUSE + R * T_STRATOSPHERE / G * Math.log(P11 / p);
const EXNER = LEVELS.map((p) => (p * 100 / P0) ** KAPPA), SPAN = EXNER.slice(1).map((e, k) => EXNER[k] - e), Z0 = LEVELS.map((p) => standard(p * 100));
const START = SPAN.map((s, k) => G * (Z0[k + 1] - Z0[k]) / (CP * s)), TOP = LEVELS.indexOf(200);
const signed = (v, unit) => `${v > 0 ? '+' : v < 0 ? '−' : ''}${Math.abs(v)} ${unit}`;
const lifts = (values) => { const out = [0]; for (let k = 0; k < K; k++) out.push(out[k] + CP * (values[k] - START[k]) * SPAN[k] / G); return out; };

export function mountExner(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), B = termColor('b');
  const theta = Float64Array.from(START), shown = Float64Array.from(START);
  let chosen = 3, scale = SCALES[0], grab = null, geometry = null;

  const fig = new Figure(root, { height: 500, minHeight: 460, step, draw });
  caption(root, 'In this plot every layer is a straight line, and its steepness is its θ: for every 0.1 that Π falls, a layer climbs about 1 km for each 100 K of its θ. The blue numbers are each layer’s θ, starting from the standard atmosphere. The lower panel shows on a finer scale how far each surface has moved, and there each layer’s steepness is how much you warmed or cooled it.');
  legend(root, [['line', 'the height of each pressure surface, and how far it has moved', A], ['faint', 'the standard atmosphere you started from', LINE]]);
  const layer = slider(controls, { label: 'Layer', min: 1, max: K, step: 1, value: chosen + 1, format: (v) => `${LEVELS[v - 1]} to ${LEVELS[v]} hPa`, onInput: (v) => choose(v - 1) });
  const warm = slider(controls, { label: 'Warm or cool it by', min: -RANGE, max: RANGE, step: 1, value: 0, format: (v) => signed(v, 'K'), onInput: set });
  buttons(controls, [['Reset', () => { theta.set(START); warm.value = 0; fig.play(true); }]]);
  const out = readout(controls);

  function choose(k) { chosen = k; layer.value = k + 1; warm.value = Math.round(theta[k] - START[k]); fig.render(); }
  function set(v) { const d = clamp(Math.round(v), -RANGE, RANGE); theta[chosen] = shown[chosen] = START[chosen] + d; warm.value = d; fig.render(); fig.play(true); }

  function step(dt) {
    const a = 1 - Math.exp(-dt * 10), most = Math.max(...lifts(theta).map(Math.abs)), target = SCALES.find((s) => most <= 0.95 * s) ?? SCALES.at(-1);
    let moving = false, changed = false;
    const ease = (v, to, eps) => { const next = Math.abs(to - v) > eps ? v + (to - v) * a : to; moving ||= next !== to; changed ||= next !== v; return next; };
    for (let k = 0; k < K; k++) shown[k] = ease(shown[k], theta[k], 0.02);
    scale = ease(scale, target, 0.5);
    if (!moving) fig.play(false);
    return changed;
  }

  function draw(ctx, w, h) {
    const left = 44, right = w - 16, top = 46, bottom = h - 38, k = chosen, ground = bottom - 30 - Math.round((bottom - top - 30) * 0.3), lower = ground + 30, mid = (lower + bottom) / 2;
    const x = (e) => left + (1 - e) / (1 - EDGE) * (right - left), y = (m) => ground - m / ZMAX * (ground - top), dy = (m) => mid - m / scale * (bottom - lower) / 2;
    const lift = lifts(shown), z = Z0.map((v, n) => v + lift[n]);
    const pts = EXNER.map((e, n) => [x(e), y(z[n])]), moved = EXNER.map((e, n) => [x(e), dy(lift[n])]);
    geometry = { lines: [pts, moved], top };
    ctx.lineWidth = 1;
    for (let km = 0; km <= ZMAX / 1000; km += 2) {
      ctx.strokeStyle = km ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, y(km * 1000)); ctx.lineTo(right, y(km * 1000)); ctx.stroke();
      if (km % 4 === 0) text(ctx, `${km} km`, left - 6, y(km * 1000), { align: 'right', color: MUTED, size: 10 });
    }
    const meters = (v) => Math.abs(v) >= 1000 ? signed(v / 1000, 'km') : signed(v, 'm'), ticks = [0, ...SCALES.flatMap((s) => [s, -s])].filter((v) => Math.abs(v) <= scale * 1.001 && (!v || Math.abs(v) >= scale / 3));
    for (const v of ticks) {
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, dy(v)); ctx.lineTo(right, dy(v)); ctx.stroke();
      text(ctx, meters(v), left - 6, dy(v), { align: 'right', color: MUTED, size: 10 });
    }
    ctx.strokeStyle = GRID; ctx.setLineDash([2, 4]);
    for (const e of EXNER) { ctx.beginPath(); ctx.moveTo(x(e), top); ctx.lineTo(x(e), ground); ctx.moveTo(x(e), lower); ctx.lineTo(x(e), bottom); ctx.stroke(); }
    ctx.setLineDash([]);
    text(ctx, 'pressure in hPa', (left + right) / 2, 14, { align: 'center', color: MUTED, size: 10 });
    LEVELS.forEach((p, n) => text(ctx, `${p}`, x(EXNER[n]), top - 10, { align: 'center', color: MUTED, size: 10 }));
    for (let e = 10; e >= EDGE * 10; e--) text(ctx, (e / 10).toFixed(1), x(e / 10), bottom + 16, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'Exner function Π, smaller higher up', (left + right) / 2, bottom + 30, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'how far each pressure surface has moved', left, ground + 17, { color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, ground, right - left, 4);

    const hint = 'drag a layer to warm or cool it', label = (n) => `${shown[n].toFixed(0)} K`, thick = `${Math.round(z[k + 1] - z[k]).toLocaleString('en-US')} m`;
    const placer = (points, segs, y0, y1, placed = []) => {
      const box = (str, size, weight, [sx, sy, f]) => {
        ctx.font = `${weight} ${size}px system-ui, sans-serif`;
        const tw = ctx.measureText(str).width, half = size * 0.6, x0 = clamp(sx - f * tw, 2, w - 2 - tw), cy = clamp(sy, y0 + half, y1 - 2 - half);
        return { x0, x1: x0 + tw, y0: cy - half, y1: cy + half };
      };
      const over = (b, o) => Math.max(0, Math.min(b.x1, o.x1) - Math.max(b.x0, o.x0)) * Math.max(0, Math.min(b.y1, o.y1) - Math.max(b.y0, o.y0));
      const cross = (b) => segs.reduce((s, [ax, ay, bx, by]) => { const m = Math.ceil(Math.hypot(bx - ax, by - ay) / 3) || 1; for (let i = 0; i <= m; i++) { const px = ax + (bx - ax) * i / m, py = ay + (by - ay) * i / m; if (px > b.x0 - 2 && px < b.x1 + 2 && py > b.y0 - 2 && py < b.y1 + 2) s += 30; } return s; }, 0);
      const crowd = (b) => cross(b) + placed.reduce((s, o) => s + over(b, o), 0) + points.reduce((s, [px, py], n) => { const r = n === k + 1 ? 9 : 5; return s + over(b, { x0: px - r, x1: px + r, y0: py - r, y1: py + r }); }, 0);
      const place = (str, size, weight, spots, optional = false) => {
        const boxes = spots.map((spot) => box(str, size, weight, spot)), b = boxes.find((c) => !crowd(c)) ?? (optional ? null : boxes.reduce((p, c) => crowd(c) < crowd(p) ? c : p));
        if (b) placed.push(b);
        return b && [b.x0, (b.y0 + b.y1) / 2];
      };
      place.reserve = (str, size, spot) => placed.push(box(str, size, 400, spot));
      return place;
    };
    const around = (points, n) => {
      const [x0, y0] = points[n], [x1, y1] = points[n + 1], len = Math.hypot(x1 - x0, y1 - y0) || 1, nx = (y0 - y1) / len, ny = (x1 - x0) / len, mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
      return [[13, 0], [13, -10], [13, 10], [24, 0], [24, -14], [24, 14], [13, -20], [13, 20], [36, 0], [-13, 0], [-13, 10], [-13, -10], [-24, 0], [48, 0]].map(([d, t]) => [mx + d * nx - t * ny, my + d * ny + t * nx, clamp(0.5 - 1.5 * Math.sign(d) * nx, 0, 1)]);
    };
    const order = [k, ...Array.from({ length: K }, (_, i) => i).filter((i) => i !== k)], path = (points) => points.slice(1).map(([px, py], n) => [...points[n], px, py]);
    const [cx0, cy0] = pts[k], [cx1, cy1] = pts[k + 1];
    const upper = placer(pts, [...path(pts), [cx0, cy0, cx0, cy1], [cx0, cy1, cx1, cy1]], top, ground), spots = [];
    upper.reserve(hint, 11, [left + 8, top + 14, 0]);
    for (let km = 0; km <= ZMAX / 1000; km += 4) upper.reserve(`${km} km`, 10, [left - 6, y(km * 1000), 1]);
    for (const n of order) spots[n] = upper(label(n), 11, n === k ? 700 : 400, around(pts, n));
    const rise = cy0 - cy1 > 24 ? upper(thick, 10, 400, [[cx0 - 7, (cy0 + cy1) / 2, 1], [cx0 + 7, (cy0 + cy1) / 2 + 8, 0], [(cx0 + cx1) / 2, cy1 - 10, 0.5]], true) : null;
    const below = placer(moved, path(moved), lower, bottom + 2), marks = [];
    for (const v of ticks) below.reserve(meters(v), 10, [left - 6, dy(v), 1]);
    for (const n of order) if (Math.abs(shown[n] - START[n]) >= 0.5) marks[n] = below(signed(Math.round(shown[n] - START[n]), 'K'), 11, n === k ? 700 : 400, around(moved, n));

    const line = (points) => {
      ctx.strokeStyle = A; ctx.lineCap = 'round';
      for (let n = 0; n < K; n++) {
        const [x0, y0] = points[n], [x1, y1] = points[n + 1];
        if (n === k) { ctx.globalAlpha = 0.25; ctx.lineWidth = 10; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke(); ctx.globalAlpha = 1; }
        ctx.lineWidth = 2.5; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      }
      ctx.fillStyle = A;
      for (const [px, py] of points) { ctx.beginPath(); ctx.arc(px, py, 3, 0, Math.PI * 2); ctx.fill(); }
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(...points[k + 1], 7, 0, Math.PI * 2); ctx.stroke();
    };
    ctx.save(); ctx.beginPath(); ctx.rect(0, top, w, ground - top); ctx.clip();
    text(ctx, hint, left + 8, top + 14, { color: MUTED, size: 11 });
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.lineJoin = 'round'; ctx.beginPath();
    EXNER.forEach((e, n) => { if (n) ctx.lineTo(x(e), y(Z0[n])); else ctx.moveTo(x(e), y(Z0[n])); });
    ctx.stroke();
    ctx.strokeStyle = MUTED; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(cx0, cy1); ctx.lineTo(cx1, cy1); ctx.stroke(); ctx.setLineDash([]);
    ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(cx0 - 4, cy1); ctx.lineTo(cx0 + 4, cy1); ctx.moveTo(cx0, cy0); ctx.lineTo(cx0, cy1); ctx.stroke();
    line(pts);
    if (rise) text(ctx, thick, ...rise, { color: INK, size: 10, halo: HALO });
    for (let n = 0; n < K; n++) text(ctx, label(n), ...spots[n], { color: B, size: 11, weight: n === k ? 700 : 400, halo: HALO });
    ctx.restore();
    ctx.save(); ctx.beginPath(); ctx.rect(0, lower - 6, w, bottom - lower + 12); ctx.clip();
    line(moved);
    marks.forEach((at, n) => { if (at) text(ctx, signed(Math.round(shown[n] - START[n]), 'K'), ...at, { color: B, size: 11, weight: n === k ? 700 : 400, halo: HALO }); });
    ctx.restore();

    const T = shown[k] * (EXNER[k] + EXNER[k + 1]) / 2 - 273.15, n = k < TOP ? TOP : K;
    out.set([
      ['layer', `${LEVELS[k]} to ${LEVELS[k + 1]} hPa`],
      ['its θ', `${shown[k].toFixed(1)} K`],
      ['steepness', `${(CP * shown[k] / G / 1e4).toFixed(2)} km for every 0.1 that Π falls`],
      ['temperature at its middle', `${T < 0 ? '−' : ''}${Math.abs(T).toFixed(1)} °C`],
      ['thickness', thick],
      [`the ${LEVELS[n]} hPa surface`, `${(z[n] / 1000).toFixed(2)} km, ${Math.abs(lift[n]) < 0.5 ? 'where it started' : `${Math.abs(lift[n]).toFixed(0)} m ${lift[n] > 0 ? 'higher' : 'lower'} than at the start`}`],
    ]);
  }

  function pick({ x: px, y: py }) {
    if (!geometry || py < geometry.top - 4) return -1;
    let best = -1, gap = 22;
    for (const pts of geometry.lines) for (let k = 0; k < K; k++) { const d = Math.hypot(px - pts[k + 1][0], py - pts[k + 1][1]); if (d < gap) { gap = d; best = k; } }
    if (best >= 0) return best;
    gap = 26;
    for (const pts of geometry.lines) for (let k = 0; k < K; k++) {
      const [x0, y0] = pts[k], [x1, y1] = pts[k + 1];
      if (px < x0 - 2 || px > x1 + 2) continue;
      const d = Math.abs(py - (y0 + (y1 - y0) * clamp((px - x0) / (x1 - x0), 0, 1)));
      if (d < gap) { gap = d; best = k; }
    }
    return best;
  }

  fig.pointer({
    hit: (p) => pick(p) >= 0,
    down: (p) => { choose(pick(p)); grab = { y: p.y, from: theta[chosen] - START[chosen] }; },
    move: (p) => { if (grab) set(grab.from + (grab.y - p.y) * GAIN); },
    up: () => { grab = null; },
  });
  fig.canvas.addEventListener('pointermove', (e) => { if (grab) return; const r = fig.canvas.getBoundingClientRect(); fig.canvas.style.cursor = pick({ x: e.clientX - r.left, y: e.clientY - r.top }) >= 0 ? 'ns-resize' : ''; });
}
