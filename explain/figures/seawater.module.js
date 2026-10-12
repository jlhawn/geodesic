import { Figure, choice, readout, legend, caption, text, rampRGB, clamp, INK, MUTED, GRID, LINE, ACCENT, CLOUD, DARK_NEUTRAL, paletteVersion } from '../runtime.module.js';
import { salinityForDensity, labelTemperature } from '../../js/ocean/seawater.module.js';
import { classPoints, water, nearestClass, mixture, classesByLatitude, columnAt, layerAt, thickClasses, nearestLatitude, ALL_CLASSES, FOUR_CLASSES, FREEZING_POINT, THIN, CELSIUS } from '../seaClasses.module.js';

const WIDE = 520, S_LO = 30, S_HI = 40, T_LO = -2, T_HI = 32, RHO_MID = 1024.8, RHO_SCALE = 1.3, LEVELS = Array.from({ length: 15 }, (_, n) => 1017 + n);
const SHALLOW = 1000, FLOOR = 4000, STRETCH = 4 / 7, HALO = 'rgba(20, 20, 22, 0.85)', DARK = `rgb(${DARK_NEUTRAL.join(', ')})`;
const ALL = classPoints(ALL_CLASSES), FOUR = classPoints(FOUR_CLASSES);
const PAIR = [{ t: CELSIUS, s: salinityForDensity(1026.8, CELSIUS) }, { t: CELSIUS + 15, s: salinityForDensity(1026.8, CELSIUS + 15) }];
const GROUPS = [['warm pool', 1021], ['thermocline', 1024.75], ['deep and bottom water', 1026.8]];

const minus = (s) => s.replace(/^-/, '−');
const celsius = (t, d = 1) => `${minus((t - CELSIUS).toFixed(d)).replace(/^−(?=0\.0*$)/, '')} °C`;
const psu = (s, d = 2) => `${s.toFixed(d)} psu`;
const density = (r, d = 2) => `${r.toFixed(d)} kg/m³`;
const signed = (v, d = 3) => `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`;
const label = (r) => (Number.isInteger(r) ? r.toFixed(1) : `${r}`);
const trim = (r) => label(Number(r.toFixed(3)));
const meters = (z) => `${Math.round(z).toLocaleString('en-US')} m`;
const degrees = (lat) => { const a = Math.abs(lat); return a < 0.05 ? 'the equator' : `${a.toFixed(1)}°${lat > 0 ? 'N' : 'S'}`; };
const shadeOf = (rho, out) => rampRGB(0.5 / (1 + Math.exp(-(rho - RHO_MID) / RHO_SCALE)), out, DARK_NEUTRAL);
const css = (c, a = 1) => `rgba(${Math.round(255 * c[0])}, ${Math.round(255 * c[1])}, ${Math.round(255 * c[2])}, ${a})`;
const isopycnal = (rho) => { const out = []; for (let t = T_LO; t <= T_HI + 1e-9; t += 0.25) out.push([salinityForDensity(rho, t + CELSIUS), t]); return out; };
const LINES = new Map(LEVELS.map((r) => [r, isopycnal(r)]));

export function mountSeawater(root) {
  const controls = root.querySelector('.controls');
  const state = { view: 'plane', probe: 'one', classes: 'all', water: { t: CELSIUS + 25, s: 35 }, waters: PAIR.map((w) => ({ ...w })), at: { lat: 0, z: 1000 } };
  let geometry = null, shading = null, sections = null, building = false, failed = false, touched = false;
  const rgb = [0, 0, 0];

  const fig = new Figure(root, { height: (w) => (w < WIDE ? Math.round(clamp(w * 1.12, 360, 420)) : 440), draw });
  const words = caption(root, '');
  const key = legend(root, []);
  choice(controls, { label: 'Show', options: [['temperature and salinity', 'plane'], ['the classes by latitude', 'section']], value: state.view, onChange: (v) => { state.view = v; arrange(); } });
  choice(controls, { label: 'Probe', options: [['one water', 'one'], ['mix two waters', 'mix']], value: state.probe, onChange: (v) => { state.probe = v; touched = false; arrange(); } });
  const probeRow = controls.lastElementChild;
  choice(controls, { label: 'Classes', options: [['all 44', 'all'], ['the four of the next figures', 'four']], value: state.classes, onChange: (v) => { state.classes = v; arrange(); } });
  const first = readout(controls), second = readout(controls), third = readout(controls);

  function arrange() {
    const plane = state.view === 'plane', four = state.classes === 'four';
    probeRow.style.display = plane ? '' : 'none';
    if (plane) {
      words.textContent = `The model’s equation of state for seawater on a plane of salinity, left to right, and temperature, up the page, shaded by density, with lines of equal density every 1 kg/m³. The dots are ${four ? 'the four classes that the next figures use, among the model’s 44' : 'the 44 classes of the model’s interior'}, each at the temperature and salinity of its label. ${state.probe === 'one' ? 'Hover or tap to move the probe, or drag it.' : 'Drag either water. The dot between them is the two mixed half and half, which averages their temperatures and salinities.'}`;
      const items = [['ramp', 'density, lighter water brighter', 'cool', DARK], ['faint', 'equal density, every 1 kg/m³', LINE], ['dash', 'the model’s freezing point', CLOUD], ['dots', four ? 'the four classes' : 'the 44 classes', four ? ACCENT : INK]];
      if (state.probe === 'mix') items.push(['dash', 'the density of each water', '#fff'], ['dash', 'the density of the mix', ACCENT]);
      else items.push(['dash', 'the density of the probe’s water; the nearest class is ringed', ACCENT]);
      key.set(items);
    } else {
      words.textContent = `The start the model builds when it has no ocean atlas to begin from, here for a sea 4000 m deep with no land, under the climate model’s opening sea surface: 28 °C on the equator and ice poleward of 72°N and 68°S. It is the same at every longitude, so each latitude has one column. Each band is ${four ? 'one of the four classes of the next figures, the same start built with those four instead of 44' : 'one of the 44 classes'}; a class thinner than 5 m is not drawn. The top kilometer is drawn four times as tall for each meter as the water below it. Hover or tap to read a class, or drag the probe.`;
      key.set([['ramp', 'density, lighter water brighter', 'cool', DARK], ['faint', 'boundaries between classes', LINE], ['line', 'the base of the mixed layer', ACCENT], ...(four ? [] : [['dash', 'the top of the 1024.0 class', '#fff']]), ['bar', 'sea ice', CLOUD]]);
      if (!sections && !building) build();
    }
    report();
    fig.render();
  }

  function build() {
    building = true;
    fig.render();
    setTimeout(() => {
      try {
        const [all, four] = classesByLatitude(16);
        sections = { all, four };
      } catch (error) { failed = true; console.info('seawater figure', error.message); }
      building = false;
      report();
      fig.render();
    }, 30);
  }

  function report() {
    if (state.view === 'plane') reportPlane(); else reportSection();
  }

  function reportPlane() {
    const points = state.classes === 'four' ? FOUR : ALL;
    if (state.probe === 'one') {
      const w = water(state.water.t, state.water.s), k = nearestClass(w.rho, points), c = points[k];
      first.set([['water', `${celsius(w.t)}, ${psu(w.s)}`], ['density', density(w.rho)], ...(w.t < FREEZING_POINT - 1e-6 ? [['below', 'the model’s freezing point']] : [])]);
      second.set([['thermal expansion', `${(w.alpha * 1e4).toFixed(2)} × 10⁻⁴ per K: ${signed(w.perKelvin)} kg/m³ for each kelvin`], ['haline contraction', `${(w.beta * 1e4).toFixed(2)} × 10⁻⁴ per psu: ${signed(w.perPsu)} kg/m³ for each psu`]]);
      const below = k > 0 ? 0.5 * (points[k - 1].rho + c.rho) : null, above = k < points.length - 1 ? 0.5 * (c.rho + points[k + 1].rho) : null;
      const span = below === null ? `water lighter than ${trim(above)}` : above === null ? `water denser than ${trim(below)}` : `water from ${trim(below)} to ${trim(above)}`;
      third.set([['nearest class', `${label(c.rho)} kg/m³, labeled ${celsius(c.t, 2)} at ${psu(c.s)}; it stands for ${span} kg/m³`]]);
    } else {
      const [a, b] = state.waters, m = mixture(a, b);
      const ra = water(a.t, a.s).rho, rb = water(b.t, b.s).rho;
      first.set([['water 1', `${celsius(a.t)}, ${psu(a.s)}, ${density(ra, 3)}`], ['water 2', `${celsius(b.t)}, ${psu(b.s)}, ${density(rb, 3)}`]]);
      second.set([['half and half', `${celsius(m.t, 2)}, ${psu(m.s)}, ${density(m.rho, 3)}`]]);
      third.set([['the mix against the two waters’ average', `${signed(m.excess)} kg/m³`], ['against the denser water', `${signed(m.overDenser)} kg/m³`]]);
    }
  }

  function reportSection() {
    const s = sections?.[state.classes];
    if (!s) { first.set([['the classes by latitude', failed ? 'could not be built' : 'building the start…']]); second.set([]); third.set([]); return; }
    const j = nearestLatitude(s, state.at.lat), L = s.layers, layer = layerAt(s, j, state.at.z), T0 = s.temperature[j * L], S0 = s.salinity[j * L];
    const mixed = water(T0, S0).rho;
    if (!layer) first.set([['at', `${degrees(s.latitude[j])}, ${meters(state.at.z)}, below the sea floor`]]);
    else if (layer.k === 0) first.set([['at', `${degrees(s.latitude[j])}, ${meters(state.at.z)}`], ['in', `the mixed layer, 0 to ${meters(layer.base)}`]]);
    else first.set([['at', `${degrees(s.latitude[j])}, ${meters(state.at.z)}`], ['class', `${label(s.densities[layer.k])} kg/m³, ${meters(layer.top)} to ${meters(layer.base)} down`]]);
    if (layer && layer.k > 0) {
      const n = j * L + layer.k, points = state.classes === 'four' ? FOUR : ALL, c = points[layer.k - 1];
      second.set([['its water here', `${celsius(s.temperature[n], 2)}, ${psu(s.salinity[n])}`], ['its label', `${celsius(c.t, 2)}, ${psu(c.s)}`]]);
    } else second.set([['the mixed layer here', `${celsius(T0, 2)}, ${psu(S0)}, ${density(mixed)}`]]);
    const thick = thickClasses(s, j);
    third.set([['this column', `the mixed layer over ${thick} ${thick === 1 ? 'class' : 'classes'} thicker than 5 m`], ['sea surface', `${signed(s.eta[j], 2)} m`]]);
  }

  function layout(w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11, top = 32, bottom = h - 38;
    return state.view === 'plane'
      ? { narrow, size, left: narrow ? 34 : 44, right: w - (narrow ? 12 : 16), top, bottom }
      : { narrow, size, left: narrow ? 46 : 54, right: w - (narrow ? 10 : 16), top, bottom };
  }

  function draw(ctx, w, h) {
    const g = layout(w, h);
    if (state.view === 'plane') plane(ctx, w, h, g); else section(ctx, w, h, g);
  }

  function paintShading(g) {
    const pw = Math.round(g.right - g.left), ph = Math.round(g.bottom - g.top), keyOf = `${pw}x${ph}:${paletteVersion}`;
    if (shading?.key === keyOf) return shading.canvas;
    const cols = Math.ceil(pw / 2), rows = Math.ceil(ph / 2), canvas = document.createElement('canvas');
    canvas.width = cols; canvas.height = rows;
    const c2 = canvas.getContext('2d'), image = c2.createImageData(cols, rows);
    for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
      const s = S_LO + (c + 0.5) / cols * (S_HI - S_LO), t = T_HI - (r + 0.5) / rows * (T_HI - T_LO);
      shadeOf(water(t + CELSIUS, s).rho, rgb);
      const n = 4 * (r * cols + c);
      image.data[n] = 255 * rgb[0]; image.data[n + 1] = 255 * rgb[1]; image.data[n + 2] = 255 * rgb[2]; image.data[n + 3] = 235;
    }
    c2.putImageData(image, 0, 0);
    shading = { key: keyOf, canvas };
    return canvas;
  }

  function plane(ctx, w, h, g) {
    const { narrow, size, left, right, top, bottom } = g, four = state.classes === 'four';
    const x = (s) => left + (s - S_LO) / (S_HI - S_LO) * (right - left), y = (t) => bottom - (t - T_LO) / (T_HI - T_LO) * (bottom - top);
    const sOf = (px) => S_LO + (px - left) / (right - left) * (S_HI - S_LO), tOf = (py) => T_LO + (bottom - py) / (bottom - top) * (T_HI - T_LO);
    geometry = { kind: 'plane', left, right, top, bottom, x, y, sOf, tOf };
    text(ctx, narrow ? 'the density of seawater' : 'the density of seawater, in the model’s equation of state', left, top - 18, { color: INK, size: size + 0.5, weight: 500 });
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(paintShading(g), left, top, right - left, bottom - top);

    ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    for (const line of LINES.values()) curve(ctx, line, x, y);
    ctx.strokeStyle = CLOUD; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(left, y(FREEZING_POINT - CELSIUS)); ctx.lineTo(right, y(FREEZING_POINT - CELSIUS)); ctx.stroke(); ctx.setLineDash([]);
    ctx.restore();

    for (const [r] of LINES) {
      if (narrow && r % 2) continue;
      const sTop = salinityForDensity(r, T_HI + CELSIUS);
      if (sTop >= S_LO && sTop <= S_HI - 0.3) text(ctx, `${r}`, x(sTop) + 2, top + 8, { color: '#fff', size: size - 0.5, halo: HALO });
      else if (sTop > S_HI) { const t = labelTemperature(r, S_HI) - CELSIUS; if (t > T_LO + 1.5 && t < T_HI - 1) text(ctx, `${r}`, right - 3, y(t) + 8, { align: 'right', color: '#fff', size: size - 0.5, halo: HALO }); }
    }
    text(ctx, narrow ? 'freezing, −1.8 °C' : 'the model’s freezing point, −1.8 °C', left + 5, y(FREEZING_POINT - CELSIUS) - 8, { color: CLOUD, size, halo: HALO });

    ctx.lineWidth = 1;
    for (let t = 0; t <= 30; t += narrow ? 10 : 5) { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left - 4, y(t)); ctx.lineTo(left, y(t)); ctx.stroke(); text(ctx, `${t}`, left - 6, y(t), { align: 'right', color: MUTED, size }); }
    text(ctx, '°C', left - 6, top - 4, { align: 'right', color: MUTED, size, baseline: 'bottom' });
    for (let s = S_LO; s <= S_HI; s += narrow ? 2 : 1) { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x(s), bottom); ctx.lineTo(x(s), bottom + 4); ctx.stroke(); text(ctx, `${s}`, x(s), bottom + 13, { align: 'center', color: MUTED, size }); }
    text(ctx, 'salinity, psu', (left + right) / 2, bottom + 28, { align: 'center', color: MUTED, size });

    const points = four ? FOUR : ALL;
    for (const p of ALL) {
      const lit = !four || FOUR.some((q) => q.rho === p.rho);
      if (lit && four) continue;
      dot(ctx, x(p.s), y(p.t - CELSIUS), narrow ? 1.9 : 2.3, four ? 'rgba(221, 221, 221, 0.3)' : INK);
    }
    if (four) for (const p of FOUR) { dot(ctx, x(p.s), y(p.t - CELSIUS), 4, ACCENT); text(ctx, label(p.rho), x(p.s) + 8, y(p.t - CELSIUS), { color: ACCENT, size, halo: HALO }); }
    else if (state.probe === 'one') groups(ctx, g, x, y);

    if (state.probe === 'one') {
      const w0 = water(state.water.t, state.water.s), k = nearestClass(w0.rho, points), c = points[k];
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]); curve(ctx, isopycnal(w0.rho), x, y); ctx.setLineDash([]);
      ctx.restore();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(x(c.s), y(c.t - CELSIUS), 5.5, 0, Math.PI * 2); ctx.stroke();
      const px = x(w0.s), py = y(w0.t - CELSIUS);
      ring(ctx, px, py);
      const right_ = px > (left + right) / 2, words2 = density(w0.rho);
      text(ctx, words2, right_ ? px - 12 : px + 12, py - 12, { align: right_ ? 'right' : 'left', color: '#fff', size, weight: 500, halo: HALO });
      geometry.handles = [[px, py]];
      if (!touched) {
        ctx.font = `400 ${size}px system-ui, sans-serif`;
        const taken = four ? FOUR.map((p) => { const lx = x(p.s) + 8, ly = y(p.t - CELSIUS); return [lx - 4, ly - 7, lx + ctx.measureText(label(p.rho)).width + 2, ly + 7]; }) : [];
        ctx.font = '400 11px system-ui, sans-serif';
        const width = ctx.measureText('drag').width, spots = [[right_, py + 12], [!right_, py + 12], [!right_, py - 12]];
        const box = ([leftward, hy]) => (leftward ? [px - 12 - width, hy - 7, px - 12, hy + 7] : [px + 12, hy - 7, px + 12 + width, hy + 7]);
        const free = spots.find((spot) => { const [a0, b0, a1, b1] = box(spot); return taken.every(([c0, d0, c1, d1]) => a1 < c0 || a0 > c1 || b1 < d0 || b0 > d1); }) ?? spots[0];
        text(ctx, 'drag', free[0] ? px - 12 : px + 12, free[1], { align: free[0] ? 'right' : 'left', color: ACCENT, size: 11, halo: HALO });
      }
    } else {
      const [a, b] = state.waters, m = mixture(a, b);
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
      ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.85)';
      for (const v of state.waters) curve(ctx, isopycnal(water(v.t, v.s).rho), x, y);
      ctx.strokeStyle = ACCENT; curve(ctx, isopycnal(m.rho), x, y);
      ctx.setLineDash([]);
      ctx.restore();
      ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(x(a.s), y(a.t - CELSIUS)); ctx.lineTo(x(b.s), y(b.t - CELSIUS)); ctx.stroke();
      dot(ctx, x(m.s), y(m.t - CELSIUS), 4.5, ACCENT);
      const mx = x(m.s), my = y(m.t - CELSIUS);
      ctx.font = `500 ${size}px system-ui, sans-serif`;
      const side = mx + 10 + ctx.measureText('half and half').width > right - 4;
      text(ctx, 'half and half', side ? mx - 10 : mx + 10, my + 1, { align: side ? 'right' : 'left', color: ACCENT, size, weight: 500, halo: HALO });
      geometry.handles = state.waters.map((v, n) => {
        const px = x(v.s), py = y(v.t - CELSIUS);
        ring(ctx, px, py);
        text(ctx, `${n + 1}`, px, py - 14, { align: 'center', color: '#fff', size: size + 1, weight: 600, halo: HALO });
        return [px, py];
      });
      if (!touched) { const [px, py] = geometry.handles[1]; text(ctx, 'drag', px + 12, py + 1, { color: ACCENT, size: 11, halo: HALO }); }
    }
  }

  function groups(ctx, g, x, y) {
    const { narrow, size, left, right, top, bottom } = g;
    for (const [name, rho] of GROUPS) {
      if (narrow && name === 'warm pool') continue;
      const p = ALL.find((q) => q.rho === rho), px = x(p.s), py = y(p.t - CELSIUS);
      const words = narrow && name === 'deep and bottom water' ? 'deep water' : name;
      ctx.font = `400 ${size}px system-ui, sans-serif`;
      const width = ctx.measureText(words).width, leftward = px + 10 + width > right - 4;
      const ly = clamp(name === 'deep and bottom water' ? py + 2 : py, top + 10, bottom - 10);
      text(ctx, words, leftward ? px - 10 : px + 10, ly, { align: leftward ? 'right' : 'left', color: INK, size, halo: HALO });
    }
    const outflow = ALL.slice(-3), words = narrow ? 'Red Sea, Mediterranean' : 'Red Sea and Mediterranean';
    ctx.font = `400 ${size}px system-ui, sans-serif`;
    const width = ctx.measureText(words).width, cx = clamp(outflow.reduce((sum, p) => sum + x(p.s), 0) / 3, left + width / 2 + 3, right - width / 2 - (narrow ? 3 : 34));
    text(ctx, words, cx, y(Math.min(...outflow.map((p) => p.t - CELSIUS))) + 14, { align: 'center', color: INK, size, halo: HALO });
  }

  function section(ctx, w, h, g) {
    const { narrow, size, left, right, top, bottom } = g, s = sections?.[state.classes];
    const x = (lat) => left + (lat + 90) / 180 * (right - left), H1 = STRETCH * (bottom - top), H2 = bottom - top - H1;
    const y = (z) => (z <= SHALLOW ? top + z / SHALLOW * H1 : top + H1 + (Math.min(z, FLOOR) - SHALLOW) / (FLOOR - SHALLOW) * H2);
    const zOf = (py) => { const d = py - top; return d <= H1 ? d / H1 * SHALLOW : SHALLOW + (d - H1) / H2 * (FLOOR - SHALLOW); };
    geometry = { kind: 'section', left, right, top, bottom, x, y, latOf: (px) => -90 + (px - left) / (right - left) * 180, zOf };
    text(ctx, narrow ? 'the classes at the start' : 'the classes at the start, by latitude', left, top - 18, { color: INK, size: size + 0.5, weight: 500 });
    if (!s) { text(ctx, failed ? 'The start could not be built.' : 'Building the start with the model’s ocean…', (left + right) / 2, (top + bottom) / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    const L = s.layers, n = s.latitude.length, bases = Array.from({ length: n }, (_, j) => columnAt(s, j));
    const highlight = (() => { const j = nearestLatitude(s, state.at.lat), layer = layerAt(s, j, state.at.z); return layer ? layer.k : -1; })();

    for (let k = 1; k < L; k++) {
      shadeOf(s.densities[k], rgb);
      ctx.fillStyle = css(rgb, 0.92);
      band(ctx, bases, k, x, y, s);
      ctx.fill();
    }
    ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(255, 255, 255, 0.16)';
    for (let k = 1; k < L - 1; k++) edge(ctx, s, bases, k, x, y, (j) => (s.thickness[j * L + k] > THIN || s.thickness[j * L + k + 1] > THIN) && bases[j][k] < FLOOR - 1);
    for (let j = 0; j < n; j++) {
      const a = x(j ? 0.5 * (s.latitude[j - 1] + s.latitude[j]) : -90), b = x(j < n - 1 ? 0.5 * (s.latitude[j] + s.latitude[j + 1]) : 90);
      shadeOf(water(s.temperature[j * L], s.salinity[j * L]).rho, rgb);
      ctx.fillStyle = css(rgb, 0.92); ctx.fillRect(a, top, b - a + 0.5, y(bases[j][0]) - top);
    }
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath();
    for (let j = 0; j < n; j++) { const px = x(s.latitude[j]), py = y(bases[j][0]); if (j) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    ctx.stroke();
    ctx.fillStyle = CLOUD;
    for (let j = 0; j < n; j++) if (s.iced[j]) { const a = x(j ? 0.5 * (s.latitude[j - 1] + s.latitude[j]) : -90), b = x(j < n - 1 ? 0.5 * (s.latitude[j] + s.latitude[j + 1]) : 90); ctx.fillRect(a, top - 4, b - a + 0.5, 4); }

    const marker = state.classes === 'four' ? -1 : s.densities.indexOf(1024);
    if (marker > 0) {
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]); ctx.beginPath();
      let on = false, peak = null;
      for (let j = 0; j < n; j++) {
        const thick = s.thickness[j * L + marker] > THIN, px = x(s.latitude[j]), py = y(bases[j][marker - 1]);
        if (thick) { if (on) ctx.lineTo(px, py); else ctx.moveTo(px, py); on = true; if (!peak || Math.abs(s.latitude[j]) < Math.abs(peak[2])) peak = [px, py, s.latitude[j], bases[j][marker - 1]]; } else on = false;
      }
      ctx.stroke(); ctx.setLineDash([]);
      if (peak) text(ctx, `${narrow ? '' : 'top of '}1024.0: ${meters(peak[3])}`, peak[0] + 6, peak[1] - 9, { color: '#fff', size, halo: HALO });
    }
    if (highlight > 0) {
      ctx.fillStyle = 'rgba(255, 232, 160, 0.12)';
      band(ctx, bases, highlight, x, y, s); ctx.fill();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5;
      const thick = (j) => s.thickness[j * L + highlight] > THIN;
      edge(ctx, s, bases, highlight - 1, x, y, thick);
      edge(ctx, s, bases, highlight, x, y, thick);
    }
    text(ctx, narrow ? 'mixed layer' : 'the mixed layer, 60 m', x(narrow ? -86 : -86), y(30) + 0.5, { color: '#fff', size, halo: HALO });

    ctx.lineWidth = 1;
    const depths = narrow ? [0, 250, 500, 750, 1000, 2000, 3000, 4000] : [0, 200, 400, 600, 800, 1000, 2000, 3000, 4000];
    for (const z of depths) { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left - 4, y(z)); ctx.lineTo(left, y(z)); ctx.stroke(); text(ctx, z ? z.toLocaleString('en-US') : '0', left - 6, y(z), { align: 'right', color: MUTED, size }); }
    text(ctx, 'm', left - 6, top - 18, { align: 'right', color: MUTED, size });
    ctx.strokeStyle = GRID; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(left - 10, y(SHALLOW)); ctx.lineTo(left, y(SHALLOW)); ctx.stroke(); ctx.setLineDash([]);
    for (let lat = -90; lat <= 90; lat += 30) {
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x(lat), bottom); ctx.lineTo(x(lat), bottom + 4); ctx.stroke();
      if (narrow && Math.abs(lat) === 90) continue;
      text(ctx, lat === 0 ? '0°' : `${Math.abs(lat)}°${lat > 0 ? 'N' : 'S'}`, x(lat), bottom + 13, { align: 'center', color: MUTED, size });
    }
    text(ctx, 'latitude', (left + right) / 2, bottom + 28, { align: 'center', color: MUTED, size });

    const px = x(state.at.lat), py = y(state.at.z);
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.45)'; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.stroke(); ctx.setLineDash([]);
    ring(ctx, px, py);
    geometry.handles = [[px, py]];
    if (!touched) { const side = px > (left + right) / 2; text(ctx, 'drag', side ? px - 12 : px + 12, py - 12, { align: side ? 'right' : 'left', color: ACCENT, size: 11, halo: HALO }); }
  }

  function edge(ctx, s, bases, k, x, y, keep) {
    ctx.beginPath();
    let on = false;
    for (let j = 0; j < bases.length; j++) {
      const px = x(s.latitude[j]), py = y(bases[j][k]);
      if (keep(j)) { if (on) ctx.lineTo(px, py); else ctx.moveTo(px, py); on = true; } else on = false;
    }
    ctx.stroke();
  }

  function band(ctx, bases, k, x, y, s) {
    const n = bases.length;
    ctx.beginPath();
    for (let j = 0; j < n; j++) { const px = x(s.latitude[j]), py = y(k ? bases[j][k - 1] : 0); if (j) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    for (let j = n - 1; j >= 0; j--) ctx.lineTo(x(s.latitude[j]), y(bases[j][k]));
    ctx.closePath();
  }

  function pick(at) {
    const gm = geometry;
    if (!gm || at.x < gm.left - 10 || at.x > gm.right + 10 || at.y < gm.top - 10 || at.y > gm.bottom + 10) return false;
    const px = clamp(at.x, gm.left, gm.right), py = clamp(at.y, gm.top, gm.bottom);
    touched = true;
    if (gm.kind === 'plane') {
      const s = Math.round(gm.sOf(px) * 100) / 100, t = Math.round(gm.tOf(py) * 10) / 10 + CELSIUS;
      if (state.probe === 'one') state.water = { t, s };
      else {
        const [a, b] = gm.handles, n = Math.hypot(at.x - a[0], at.y - a[1]) <= Math.hypot(at.x - b[0], at.y - b[1]) ? 0 : 1;
        state.waters[grabbed ?? n] = { t, s };
      }
    } else state.at = { lat: Math.round(gm.latOf(px) * 10) / 10, z: Math.round(gm.zOf(py)) };
    report();
    fig.render();
    return true;
  }

  let grabbed = null;
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const near = (p) => {
    const handles = geometry?.handles ?? [];
    let best = -1, d = 24;
    handles.forEach(([hx, hy], n) => { const e = Math.hypot(p.x - hx, p.y - hy); if (e < d) { d = e; best = n; } });
    return best;
  };
  fig.pointer({ hit: (p) => near(p) >= 0, down: (p) => { grabbed = state.view === 'plane' && state.probe === 'mix' ? near(p) : null; return pick(p); }, move: pick, up: () => { grabbed = null; } });
  fig.canvas.addEventListener('click', (e) => { if (grabbed === null) pick(at(e)); });
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType === 'mouse' && !e.buttons && !(state.view === 'plane' && state.probe === 'mix')) pick(at(e)); });

  arrange();
}

function curve(ctx, points, x, y) {
  ctx.beginPath();
  let on = false;
  for (const [s, t] of points) {
    if (!Number.isFinite(s)) { on = false; continue; }
    if (on) ctx.lineTo(x(s), y(t)); else ctx.moveTo(x(s), y(t));
    on = true;
  }
  ctx.stroke();
}

function dot(ctx, px, py, r, color) {
  ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 1.2;
  ctx.beginPath(); ctx.arc(px, py, r, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
}

function ring(ctx, px, py) {
  ctx.fillStyle = HALO; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(px, py, 6.5, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  ctx.fillStyle = ACCENT; ctx.beginPath(); ctx.arc(px, py, 2.5, 0, Math.PI * 2); ctx.fill();
}
