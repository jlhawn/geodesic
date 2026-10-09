import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, anomalyColor, rampColor, clamp, ACCENT, INK, MUTED, LINE, GRID, P0 } from '../runtime.module.js';
import { createSlice } from '../sliceCore.module.js';

export function mountCells(root) {
  const controls = root.querySelector('.controls');
  const slice = createSlice({ dt: 10 });
  const { M, K, heated: c, heatedHalf, dSigma } = slice;
  const z0 = Float64Array.from({ length: K + 1 }, (_, k) => slice.interfaceHeight(k, 0));
  const T0k = Float64Array.from({ length: K }, (_, k) => slice.temperature(k, 0));
  const SPEED = 600, zMax = 18000, LAPSE_LIMIT = 2e-3, EXAGGERATION = { isobars: 20, isotherms: 5 };
  const lapseAt = (j, i) => (slice.temperature(j, i) - slice.temperature(j - 1, i)) / (slice.layerHeight(j - 1, i) - slice.layerHeight(j, i));
  function tropopause(i) {
    for (let j = K - 1; j >= 1; j--) {
      if (lapseAt(j, i) >= LAPSE_LIMIT) continue;
      if (j === K - 1) return slice.interfaceHeight(j, i);
      const below = slice.interfaceHeight(j + 1, i), above = slice.interfaceHeight(j, i), lb = lapseAt(j + 1, i), la = lapseAt(j, i);
      return below + (lb - LAPSE_LIMIT) / (lb - la) * (above - below);
    }
    return null;
  }
  const trop0 = tropopause(0);
  const ISOTHERMS = [10, 0, -10, -20, -30, -40, -50];
  const iso0 = {};
  function isotherm(celsius, i) {
    const target = celsius + 273.15;
    for (let k = K - 1; k >= 1; k--) {
      const lower = slice.temperature(k, i), upper = slice.temperature(k - 1, i);
      if ((target - lower) * (target - upper) > 0) continue;
      const zl = slice.layerHeight(k, i), zu = slice.layerHeight(k - 1, i);
      return lower === upper ? zl : zl + (target - lower) / (upper - lower) * (zu - zl);
    }
    return null;
  }
  for (const celsius of ISOTHERMS) iso0[celsius] = isotherm(celsius, 0);
  let heating = 1, mode = 'column', friction = false, colors = 'delta', lines = 'isobars';
  const FILLS = {
    delta: { key: ['gradient', 'cooler or warmer than at the start', 'cool', 'warm'], color: (k, i) => anomalyColor((slice.temperature(k, i) - T0k[k]) / 2.5, 0.6) },
    temperature: { key: ['ramp', 'temperature, −60 to 30 °C', 'cool', 'warm', 'neutral'], color: (k, i) => rampColor((slice.temperature(k, i) - 213.15) / 90) },
    theta: { key: ['ramp', 'θ, 285 to 345 K', 'cool', 'warm', 'neutral'], color: (k, i) => rampColor((slice.theta[k * M + i] - 285) / 60) },
  };
  const legendFor = () => [FILLS[colors].key, ['arrow', 'wind between columns', 'rgba(255,255,255,0.85)'], ['varrow', 'rising or sinking air', ACCENT], ['bar', 'change in surface pressure', 'cool'], ['dots', 'tropopause', 'rgba(255,255,255,0.7)'], ...(lines === 'isotherms' ? [['faint', 'isotherms every 10 °C', 'rgba(255,255,255,0.6)']] : [])];
  const captionFor = () => `Height changes are drawn ${lines === 'isobars' ? 'twenty' : 'five'} times larger than they are.`;

  const fig = new Figure(root, { height: 470, minHeight: 360, step, draw });
  const note = caption(root, captionFor());
  const key = legend(root, legendFor());
  controls.classList.add('two');
  choice(controls, { label: 'Color shows', options: [['change since the start', 'delta'], ['temperature', 'temperature'], ['θ', 'theta']], value: 'delta', onChange: (v) => { colors = v; key.set(legendFor()); fig.render(); }, span: true });
  choice(controls, { label: 'Lines show', options: [['pressure surfaces', 'isobars'], ['isotherms', 'isotherms']], value: 'isobars', onChange: (v) => { lines = v; note.textContent = captionFor(); key.set(legendFor()); fig.render(); }, span: true });
  slider(controls, { label: 'Warm the middle three columns by', min: 0, max: 3, step: 0.1, value: 1, format: (v) => `${v.toFixed(1)} °C per hour`, onInput: (v) => { heating = v; }, span: true });
  choice(controls, { label: 'Heat goes into', options: [['the whole column', 'column'], ['the ground', 'ground']], value: 'column', onChange: (v) => { mode = v; }, span: true });
  choice(controls, { label: 'Friction', options: [['off', 'off'], ['on', 'on']], value: 'off', onChange: (v) => { friction = v === 'on'; } });
  buttons(controls, [['Restart', () => { slice.reset(); show(); fig.render(); }], ['Pause', (b) => { fig.play(!fig.running); b.textContent = fig.running ? 'Pause' : 'Run'; }]]);
  const out = readout(controls);

  function step(dt) {
    let remaining = dt * SPEED;
    for (let n = 0; remaining > 0 && n < 200; n++) { const h = Math.min(slice.dt, remaining); slice.step({ heating: heating / 3600, mode, friction }, h); remaining -= h; }
    show(30);
  }

  function show(every = 1) {
    let maxU = 0, maxW = 0;
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) { maxU = Math.max(maxU, Math.abs(slice.u[k * M + i])); maxW = Math.max(maxW, Math.abs(slice.verticalVelocity(k, i))); }
    const minutes = Math.round(slice.state.time / 60);
    out.set([['time', `${Math.floor(minutes / 60)} h ${String(minutes % 60).padStart(2, '0')} min`], ['strongest wind', `${maxU.toFixed(1)} m/s`], ['strongest updraught', `${(maxW * 100).toFixed(0)} cm/s`]], every);
  }

  function draw(ctx, w, h) {
    const left = 44, right = 56, top = 18, barsH = 90, bottom = h - barsH - 24;
    const visible = w >= 560 ? 11 : 7, half = (visible - 1) / 2;
    const plotW = w - left - right, colW = plotW / visible;
    const xc = (i) => left + (i - (c - half) + 0.5) * colW;
    const yz = (z) => bottom - (z / zMax) * (bottom - top);
    const EXAG = EXAGGERATION[lines];
    const zi = (k, i) => z0[k] + EXAG * (slice.interfaceHeight(k, i) - z0[k]);
    const yEdge = (k, i, side) => yz(0.5 * (zi(k, i) + zi(k, i + side)));
    const yTop = (k, i, side) => k === 0 ? top - 5 : side === 0 ? yz(zi(k, i)) : yEdge(k, i, side);
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, plotW, bottom - top); ctx.clip();
    for (let i = c - half - 1; i <= c + half + 1; i++) {
      const xl = xc(i) - colW / 2, x = xc(i), xr = xc(i) + colW / 2;
      for (let k = 0; k < K; k++) {
        ctx.fillStyle = FILLS[colors].color(k, i);
        ctx.globalAlpha = colors === 'delta' ? 1 : 0.6;
        ctx.beginPath();
        ctx.moveTo(xl, yTop(k, i, -1)); ctx.lineTo(x, yTop(k, i, 0)); ctx.lineTo(xr, yTop(k, i, 1));
        ctx.lineTo(xr, yEdge(k + 1, i, 1)); ctx.lineTo(x, yz(zi(k + 1, i))); ctx.lineTo(xl, yEdge(k + 1, i, -1));
        ctx.closePath(); ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    if (lines === 'isobars') for (let k = 1; k <= K; k++) { ctx.beginPath(); for (let i = c - half - 1; i <= c + half + 1; i++) { const x = xc(i), y = yz(zi(k, i)); if (i === c - half - 1) ctx.moveTo(x, y); else ctx.lineTo(x, y); } ctx.stroke(); }
    ctx.strokeStyle = GRID;
    for (let i = c - half; i <= c + half + 1; i++) { const x = xc(i) - colW / 2; ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, bottom); ctx.stroke(); }
    const uScale = colW * 0.5 / 3;
    for (let e = c - half - 1; e <= c + half; e++) {
      const x = xc(e) + colW / 2;
      for (let k = 1; k < K; k++) {
        const u = slice.u[k * M + e];
        if (Math.abs(u) < 0.05) continue;
        const y = 0.5 * (yEdge(k, e, 1) + yEdge(k + 1, e, 1)), len = clamp(u * uScale, -colW * 0.7, colW * 0.7);
        arrow(ctx, x - len / 2, y, x + len / 2, y, { color: 'rgba(255,255,255,0.85)', width: 1.5, head: 5 });
      }
    }
    for (let i = c - half - 1; i <= c + half + 1; i++) for (let k = 2; k < K; k++) {
      const wv = slice.verticalVelocity(k, i);
      if (Math.abs(wv) < 0.004) continue;
      const x = xc(i), y = yz(zi(k, i)), len = clamp(wv * 300, -36, 36);
      arrow(ctx, x, y + len / 2, x, y - len / 2, { color: ACCENT, width: 1.5, head: 5 });
    }
    if (lines === 'isotherms') {
      ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 1;
      for (const celsius of ISOTHERMS) {
        ctx.beginPath();
        let drawing = false, first = null;
        for (let i = c - half - 1; i <= c + half + 1; i++) {
          const z = isotherm(celsius, i);
          if (z === null) { drawing = false; continue; }
          const zd = iso0[celsius] + EXAG * (z - iso0[celsius]);
          if (first === null && i >= c - half) first = zd;
          if (drawing) ctx.lineTo(xc(i), yz(zd)); else ctx.moveTo(xc(i), yz(zd));
          drawing = true;
        }
        ctx.stroke();
        if (first !== null && colW >= 40) text(ctx, `${celsius} °C`, left + 4, yz(first) - 6, { color: 'rgba(255,255,255,0.55)', size: 9 });
      }
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.7)'; ctx.lineWidth = 2; ctx.lineCap = 'round'; ctx.setLineDash([0.5, 4]); ctx.beginPath();
    let pen = false;
    for (let i = c - half - 1; i <= c + half + 1; i++) {
      const zt = tropopause(i);
      if (zt === null) { pen = false; continue; }
      const y = yz(trop0 + EXAG * (zt - trop0));
      if (pen) ctx.lineTo(xc(i), y); else ctx.moveTo(xc(i), y);
      pen = true;
    }
    ctx.stroke(); ctx.setLineDash([]); ctx.lineWidth = 1; ctx.lineCap = 'butt';
    ctx.restore();
    for (let km = 0; km <= 15; km += 5) text(ctx, `${km} km`, left - 6, yz(km * 1000), { align: 'right', color: MUTED, size: 11 });
    for (let k = 1; k < K; k++) if (z0[k] < zMax) text(ctx, `${Math.round(k * dSigma * 1000)} hPa`, left + plotW + 6, yz(z0[k]), { color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, plotW, 5);
    if (heating > 0) { ctx.fillStyle = ACCENT; ctx.fillRect(xc(c - heatedHalf) - colW * 0.4, bottom, (2 * heatedHalf + 1) * colW - colW * 0.2, 5); text(ctx, mode === 'ground' ? 'warmed ground' : 'warmed columns', xc(c), top + 10, { align: 'center', color: ACCENT, size: 11 }); }
    const base = bottom + 5 + (h - bottom - 5) / 2, perHPa = 16;
    const narrow = colW < 44;
    ctx.save(); ctx.translate(left - 30, base); ctx.rotate(-Math.PI / 2); text(ctx, narrow ? 'pressure, hPa' : 'surface pressure', 0, 0, { align: 'center', color: MUTED, size: 10 }); ctx.restore();
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, base); ctx.lineTo(left + plotW, base); ctx.stroke();
    for (let i = c - half; i <= c + half; i++) {
      const d = (slice.pi[i] - P0) / 100, len = clamp(d * perHPa, -36, 36);
      ctx.fillStyle = anomalyColor(Math.sign(d) * 0.8, 0.8);
      ctx.fillRect(xc(i) - colW * 0.3, Math.min(base, base - len), colW * 0.6, Math.abs(len));
      if (!narrow || i === c || Math.abs(i - c) === half) text(ctx, `${d >= 0 ? '+' : '−'}${Math.abs(d).toFixed(1)}${colW < 70 ? '' : ' hPa'}`, xc(i), d >= 0 ? base - len - 9 : base - len + 9, { align: 'center', color: INK, size: 10 });
    }
  }

  show();
  fig.play(true);
}
