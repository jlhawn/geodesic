import { Figure, choice, readout, legend, caption, text, rampRGB, clamp, INK, MUTED, LINE, ACCENT, CLOUD, DARK_NEUTRAL, paletteVersion } from '../runtime.module.js';
import { salinityForDensity, labelTemperature } from '../../js/ocean/seawater.module.js';
import { HYBRID_LAYERS, SURFACE_LAYERS, SURFACE_DEPTH, FREEZING_POINT, CELSIUS, water, mixture, labelOf, labelRange, isopycnalLine, depthFraction, depthAt, layersByLatitude, layerAt, modeCounts, emptyEverywhere, nearestLatitude } from '../seaClasses.module.js';

const WIDE = 520, S_LO = 30, S_HI = 40, T_LO = -2, T_HI = 32, RHO_MID = 1024.8, RHO_SCALE = 1.3, LEVELS = Array.from({ length: 15 }, (_, n) => 1017 + n);
const HALO = 'rgba(20, 20, 22, 0.85)', DARK = `rgb(${DARK_NEUTRAL.join(', ')})`;
const BAND = 'rgba(255, 232, 160, 0.4)', LABEL_LINE = 'rgba(255, 255, 255, 0.26)', ISOPYCNAL_EDGE = 'rgba(255, 255, 255, 0.55)', FIXED_EDGE = 'rgba(255, 255, 255, 0.34)';
const PAIR = [{ t: CELSIUS, s: salinityForDensity(1026.8, CELSIUS) }, { t: CELSIUS + 15, s: salinityForDensity(1026.8, CELSIUS + 15) }];
const DEEP = [['NADW', 2.5, 34.92], ['AABW', -0.5, 34.66]];
const DENSITIES = HYBRID_LAYERS.densities, LAYERS = HYBRID_LAYERS.K;

const minus = (s) => s.replace(/^-/, '−');
const celsius = (t, d = 1) => `${minus((t - CELSIUS).toFixed(d)).replace(/^−(?=0\.0*$)/, '')} °C`;
const psu = (s, d = 2) => `${s.toFixed(d)} psu`;
const density = (r, d = 2) => `${r.toFixed(d)} kg/m³`;
const signed = (v, d = 3) => `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`;
const label = (r) => (Number.isInteger(r) ? r.toFixed(1) : `${r}`);
const trim = (r) => label(Number(r.toFixed(3)));
const meters = (z) => `${Math.round(z).toLocaleString('en-US')} m`;
const depthText = (z) => (z < 100 ? `${z.toFixed(1)} m` : meters(z));
const degrees = (lat) => { const a = Math.abs(lat); return a < 0.05 ? 'the equator' : `${a.toFixed(1)}°${lat > 0 ? 'N' : 'S'}`; };
const ordinal = (n) => `${n}${n % 100 >= 11 && n % 100 <= 13 ? 'th' : ['th', 'st', 'nd', 'rd'][n % 10] ?? 'th'}`;
const shadeOf = (rho, out) => rampRGB(0.5 / (1 + Math.exp(-(rho - RHO_MID) / RHO_SCALE)), out, DARK_NEUTRAL);
const css = (c, a = 1) => `rgba(${Math.round(255 * c[0])}, ${Math.round(255 * c[1])}, ${Math.round(255 * c[2])}, ${a})`;
const isopycnal = (rho) => isopycnalLine(rho, T_LO, T_HI);
const LINES = new Map(LEVELS.map((r) => [r, isopycnal(r)]));
const LABEL_LINES = Array.from(DENSITIES, (r) => isopycnal(r));
const BOUNDS = Array.from({ length: LAYERS }, (_, k) => { const { below, above } = labelRange(k); return [isopycnal(Number.isFinite(below) ? below : 900), isopycnal(Number.isFinite(above) ? above : 1100)]; });

export function mountSeawater(root) {
  const controls = root.querySelector('.controls');
  const state = { view: 'plane', probe: 'one', water: { t: CELSIUS + 25, s: 35 }, waters: PAIR.map((w) => ({ ...w })), at: { lat: 0, z: 870 } };
  let geometry = null, shading = null, backdrop = null, section = null, building = false, failed = false, touched = false;
  const rgb = [0, 0, 0];

  const fig = new Figure(root, { height: (w) => (w < WIDE ? Math.round(clamp(w * 1.12, 360, 420)) : 440), draw });
  const words = caption(root, '');
  const key = legend(root, []);
  choice(controls, { label: 'Show', options: [['temperature and salinity', 'plane'], ['the layers by latitude', 'section']], value: state.view, onChange: (v) => { state.view = v; arrange(); } });
  choice(controls, { label: 'Probe', options: [['one water', 'one'], ['mix two waters', 'mix']], value: state.probe, onChange: (v) => { state.probe = v; touched = false; arrange(); } });
  const probeRow = controls.lastElementChild;
  const first = readout(controls), second = readout(controls), third = readout(controls);

  function arrange() {
    const plane = state.view === 'plane';
    probeRow.style.display = plane ? '' : 'none';
    if (plane) {
      words.textContent = `The model’s equation of state for seawater on a plane of salinity, left to right, and temperature, up the page, shaded by density, with lines of equal density every 1 kg/m³. The faint lines are the 58 label densities of the model’s layers; the twelve lightest, dotted, lie in the warm, fresh corner. NADW and AABW mark North Atlantic Deep Water and Antarctic Bottom Water. ${state.probe === 'one' ? 'Hover or tap to move the probe, or drag it. The yellow band is the water that the nearest label stands for.' : 'Drag either water. The dot between them is the two mixed half and half, which averages their temperatures and salinities.'}`;
      const items = [['ramp', 'density, lighter water brighter', 'cool', DARK], ['faint', 'equal density, every 1 kg/m³', LINE], ['faint', 'the 58 labels', LABEL_LINE], ['dots', 'the twelve surface labels', LABEL_LINE], ['dash', 'the model’s freezing point', CLOUD]];
      if (state.probe === 'mix') items.push(['dash', 'the density of each water', '#fff'], ['dash', 'the density of the mix', ACCENT]);
      else items.push(['dash', 'the density of the probe’s water', ACCENT], ['bar', 'the water of the nearest label', BAND]);
      key.set(items);
    } else {
      words.textContent = 'The start the model builds when it has no ocean atlas to begin from, here for a sea 4000 m deep with no land, under the climate model’s opening sea surface: 28 °C on the equator and ice poleward of 72°N and 68°S. It is the same at every longitude, so each latitude has one column. Each band is one of the 58 layers, shaded by the density of its water, as the first regrid leaves them. Depth is stretched toward the surface, so the twelve surface layers, 3 to 7 m each, can be seen. Hover or tap to read a layer, or drag the probe.';
      key.set([['ramp', 'density, lighter water brighter', 'cool', DARK], ['faint', 'interface of an isopycnal layer', ISOPYCNAL_EDGE], ['dots', 'interface of a fixed-depth layer', FIXED_EDGE], ['dash', 'the base of the twelve surface layers', '#fff'], ['bar', 'sea ice', CLOUD]]);
      if (!section && !building) build();
    }
    report();
    fig.render();
  }

  function build() {
    building = true;
    fig.render();
    setTimeout(() => {
      try { section = layersByLatitude(16); } catch (error) { failed = true; console.info('seawater figure', error.message); }
      building = false;
      report();
      fig.render();
    }, 30);
  }

  function report() {
    if (state.view === 'plane') reportPlane(); else reportSection();
  }

  function reportPlane() {
    if (state.probe === 'one') {
      const w = water(state.water.t, state.water.s), k = labelOf(w.rho), { below, above } = labelRange(k);
      first.set([['water', `${celsius(w.t)}, ${psu(w.s)}`], ['density', density(w.rho)], ...(w.t < FREEZING_POINT - 1e-6 ? [['below', 'the model’s freezing point']] : [])]);
      second.set([['thermal expansion', `${(w.alpha * 1e4).toFixed(2)} × 10⁻⁴ per K: ${signed(w.perKelvin)} kg/m³ for each kelvin`], ['haline contraction', `${(w.beta * 1e4).toFixed(2)} × 10⁻⁴ per psu: ${signed(w.perPsu)} kg/m³ for each psu`]]);
      const span = !Number.isFinite(below) ? `water lighter than ${trim(above)}` : !Number.isFinite(above) ? `water denser than ${trim(below)}` : `water from ${trim(below)} to ${trim(above)}`;
      const minimum = HYBRID_LAYERS.minimum[k];
      const rule = k < SURFACE_LAYERS ? 'a surface label, so its layer is fixed-depth whatever its water' : `its layer is isopycnal only where that water spans more than ${Number(minimum.toFixed(1))} m`;
      third.set([['nearest label', `${label(DENSITIES[k])} kg/m³, the ${ordinal(k + 1)} of ${LAYERS}; it stands for ${span} kg/m³; ${rule}`]]);
    } else {
      const [a, b] = state.waters, m = mixture(a, b);
      const ra = water(a.t, a.s).rho, rb = water(b.t, b.s).rho;
      first.set([['water 1', `${celsius(a.t)}, ${psu(a.s)}, ${density(ra, 3)}`], ['water 2', `${celsius(b.t)}, ${psu(b.s)}, ${density(rb, 3)}`]]);
      second.set([['half and half', `${celsius(m.t, 2)}, ${psu(m.s)}, ${density(m.rho, 3)}`]]);
      const ka = labelOf(ra), kb = labelOf(rb), km = labelOf(m.rho);
      third.set([['the mix against the two waters’ average', `${signed(m.excess)} kg/m³`], ['against the denser water', `${signed(m.overDenser)} kg/m³`], ['labels', `the waters ${ka === kb ? `share ${label(DENSITIES[ka])}` : `have ${label(DENSITIES[ka])} and ${label(DENSITIES[kb])}`}; the mix has ${label(DENSITIES[km])}`]]);
    }
  }

  function reportSection() {
    if (!section) { first.set([['the layers by latitude', failed ? 'could not be built' : 'building the start…']]); second.set([]); third.set([]); return; }
    const s = section, j = nearestLatitude(s, state.at.lat), K = s.layers, layer = layerAt(s, j, state.at.z), at = `${degrees(s.latitude[j])}, ${depthText(state.at.z)}`;
    const counts = modeCounts(s, j), column = ['this column', `${counts.fixed} fixed-depth, ${counts.isopycnal + counts.between} isopycnal and ${counts.empty} empty layers`];
    if (!layer) { first.set([['at', `${at}, below the sea floor`]]); second.set([]); third.set([column]); return; }
    const { k } = layer, n = j * K + k, rho = s.density[n], thick = layer.base - layer.top, minimum = Number(HYBRID_LAYERS.minimum[k].toFixed(1));
    first.set([['at', at], ['layer', `labeled ${label(DENSITIES[k])} kg/m³, the ${ordinal(k + 1)} of ${K}, ${depthText(layer.top)} to ${depthText(layer.base)} down`]]);
    second.set([['its water here', `${celsius(s.temperature[n], 2)}, ${psu(s.salinity[n])}, ${density(rho)}`]]);
    const gap = rho - DENSITIES[k], against = Math.abs(gap) < 0.005 ? 'at its label’s density' : `${Math.abs(gap).toFixed(2)} kg/m³ ${gap > 0 ? 'denser' : 'lighter'} than its label`;
    const mode = layer.mode === 'empty' ? 'empty: no water here has its density, so it is a sheet 0.5 m thick on the floor'
      : layer.mode === 'fixed' ? `fixed-depth: held at its minimum of ${minimum} m, with water ${against}`
      : layer.mode === 'isopycnal' ? `isopycnal: ${Number(thick.toFixed(1))} m thick against a minimum of ${minimum} m, with water ${against}`
      : `between the two: ${Number(thick.toFixed(1))} m thick against a minimum of ${minimum} m, mode weight ${layer.weight.toFixed(2)}`;
    third.set([['the layer is', mode], column]);
  }

  function layout(w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11, top = 32, bottom = h - 38;
    return state.view === 'plane'
      ? { narrow, size, left: narrow ? 34 : 44, right: w - (narrow ? 12 : 16), top, bottom }
      : { narrow, size, left: narrow ? 46 : 54, right: w - (narrow ? 10 : 16), top, bottom };
  }

  function draw(ctx, w, h) {
    const g = layout(w, h);
    if (state.view === 'plane') plane(ctx, w, h, g); else sea(ctx, w, h, g);
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
    const { narrow, size, left, right, top, bottom } = g;
    const x = (s) => left + (s - S_LO) / (S_HI - S_LO) * (right - left), y = (t) => bottom - (t - T_LO) / (T_HI - T_LO) * (bottom - top);
    const sOf = (px) => S_LO + (px - left) / (right - left) * (S_HI - S_LO), tOf = (py) => T_LO + (bottom - py) / (bottom - top) * (T_HI - T_LO);
    geometry = { kind: 'plane', left, right, top, bottom, x, y, sOf, tOf };
    text(ctx, narrow ? 'the density of seawater' : 'the density of seawater, in the model’s equation of state', left, top - 18, { color: INK, size: size + 0.5, weight: 500 });
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(paintShading(g), left, top, right - left, bottom - top);

    ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    for (const line of LINES.values()) curve(ctx, line, x, y);
    ctx.strokeStyle = LABEL_LINE; ctx.lineWidth = 0.8;
    LABEL_LINES.forEach((line, k) => { ctx.setLineDash(k < SURFACE_LAYERS ? [1.5, 3] : []); curve(ctx, line, x, y); });
    ctx.setLineDash([]);
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

    for (const [name, t, s] of DEEP) { dot(ctx, x(s), y(t), 3.2, INK); text(ctx, name, x(s) + 7, y(t) + 1, { color: INK, size, halo: HALO }); }

    if (state.probe === 'one') {
      const w0 = water(state.water.t, state.water.s), k = labelOf(w0.rho);
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
      const [lighter, denser] = BOUNDS[k];
      ctx.fillStyle = BAND; ctx.beginPath(); trace(ctx, denser, x, y, false, false); trace(ctx, lighter, x, y, true, true); ctx.closePath(); ctx.fill();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5; curve(ctx, LABEL_LINES[k], x, y);
      ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]); curve(ctx, isopycnal(w0.rho), x, y); ctx.setLineDash([]);
      ctx.restore();
      const px = x(w0.s), py = y(w0.t - CELSIUS);
      ring(ctx, px, py);
      const right_ = px > (left + right) / 2, sx = right_ ? px - 12 : px + 12, align = right_ ? 'right' : 'left';
      text(ctx, density(w0.rho), sx, py - 12, { align, color: '#fff', size, weight: 500, halo: HALO });
      text(ctx, `label ${label(DENSITIES[k])}`, sx, clamp(py - 25, top + 8, bottom), { align, color: ACCENT, size, halo: HALO });
      geometry.handles = [[px, py]];
      if (!touched) text(ctx, 'drag', sx, py + 12, { align, color: ACCENT, size: 11, halo: HALO });
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

  function paintLayers(g) {
    const pw = Math.round(g.right - g.left), ph = Math.round(g.bottom - g.top), dpr = fig.dpr || 1, keyOf = `${pw}x${ph}@${dpr}:${paletteVersion}`;
    if (backdrop?.key === keyOf && backdrop.section === section) return backdrop.canvas;
    const s = section, K = s.layers, n = s.latitude.length, canvas = document.createElement('canvas');
    canvas.width = Math.round(pw * dpr); canvas.height = Math.round(ph * dpr);
    const c2 = canvas.getContext('2d');
    c2.scale(dpr, dpr);
    const x = (lat) => (lat + 90) / 180 * pw, y = (z) => depthFraction(z) * ph;
    const edgeX = (j) => x(j ? 0.5 * (s.latitude[j - 1] + s.latitude[j]) : -90), edgeEnd = (j) => x(j < n - 1 ? 0.5 * (s.latitude[j] + s.latitude[j + 1]) : 90);
    for (let j = 0; j < n; j++) {
      const a = edgeX(j), b = edgeEnd(j), base = j * (K + 1);
      for (let k = 0; k < K; k++) {
        const z0 = s.interface[base + k], z1 = s.interface[base + k + 1];
        if (z1 - z0 < 1e-6 || z0 >= s.floor[j] - 1) continue;
        shadeOf(s.density[j * K + k], rgb);
        c2.fillStyle = css(rgb, 1);
        c2.fillRect(a, y(z0), b - a + 0.5, y(z1) - y(z0) + 0.5);
      }
    }
    c2.lineWidth = 1;
    for (const iso of [false, true]) {
      c2.strokeStyle = iso ? ISOPYCNAL_EDGE : FIXED_EDGE;
      c2.setLineDash(iso ? [] : [1.5, 3]);
      c2.beginPath();
      for (let k = 1; k < K; k++) {
        let on = false;
        for (let j = 1; j < n; j++) {
          const a = (j - 1) * (K + 1) + k, b = j * (K + 1) + k, z0 = s.interface[a], z1 = s.interface[b];
          const here = (s.weight[(j - 1) * K + k - 1] + s.weight[j * K + k - 1]) / 2 >= 0.5 === iso && !(z0 > s.floor[j - 1] - 1 && z1 > s.floor[j] - 1);
          if (!here) { on = false; continue; }
          if (!on) c2.moveTo(x(s.latitude[j - 1]), y(z0));
          c2.lineTo(x(s.latitude[j]), y(z1));
          on = true;
        }
      }
      c2.stroke();
    }
    backdrop = { key: keyOf, section, canvas };
    return canvas;
  }

  function sea(ctx, w, h, g) {
    const { narrow, size, left, right, top, bottom } = g, s = section;
    const x = (lat) => left + (lat + 90) / 180 * (right - left), y = (z) => top + depthFraction(z) * (bottom - top);
    const zOf = (py) => depthAt((py - top) / (bottom - top));
    geometry = { kind: 'section', left, right, top, bottom, x, y, latOf: (px) => -90 + (px - left) / (right - left) * 180, zOf };
    text(ctx, narrow ? 'the layers at the start' : 'the layers at the start, by latitude', left, top - 18, { color: INK, size: size + 0.5, weight: 500 });
    if (!s) { text(ctx, failed ? 'The start could not be built.' : 'Building the start with the model’s ocean…', (left + right) / 2, (top + bottom) / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    const K = s.layers, n = s.latitude.length, layer = layerAt(s, nearestLatitude(s, state.at.lat), state.at.z), highlight = layer ? layer.k : -1;
    ctx.drawImage(paintLayers(g), left, top, right - left, bottom - top);

    ctx.fillStyle = CLOUD;
    for (let j = 0; j < n; j++) if (s.iced[j]) { const a = x(j ? 0.5 * (s.latitude[j - 1] + s.latitude[j]) : -90), b = x(j < n - 1 ? 0.5 * (s.latitude[j] + s.latitude[j + 1]) : 90); ctx.fillRect(a, top - 4, b - a + 0.5, 4); }

    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.2; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(left, y(SURFACE_DEPTH)); ctx.lineTo(right, y(SURFACE_DEPTH)); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, narrow ? `surface layers, ${Math.round(SURFACE_DEPTH)} m` : `the twelve surface layers, to ${Math.round(SURFACE_DEPTH)} m`, left + 6, y(SURFACE_DEPTH) - 9, { color: '#fff', size, halo: HALO });
    const empty = emptyEverywhere(s);
    if (empty) text(ctx, narrow ? `${empty} empty layers on the floor` : `the ${empty} densest layers are empty, on the floor`, left + 6, bottom - 9, { color: '#fff', size, halo: HALO });

    if (highlight >= 0) {
      const edge = (k) => { ctx.beginPath(); for (let j = 0; j < n; j++) { const px = x(s.latitude[j]), py = y(s.interface[j * (K + 1) + k]); if (j) ctx.lineTo(px, py); else ctx.moveTo(px, py); } };
      ctx.fillStyle = 'rgba(255, 232, 160, 0.16)'; edge(highlight);
      for (let j = n - 1; j >= 0; j--) ctx.lineTo(x(s.latitude[j]), y(s.interface[j * (K + 1) + highlight + 1]));
      ctx.closePath(); ctx.fill();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5;
      edge(highlight); ctx.stroke();
      edge(highlight + 1); ctx.stroke();
    }

    ctx.lineWidth = 1;
    const depths = narrow ? [0, 20, 50, 100, 200, 500, 1000, 2000, 4000] : [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000];
    for (const z of depths) { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left - 4, y(z)); ctx.lineTo(left, y(z)); ctx.stroke(); text(ctx, z ? z.toLocaleString('en-US') : '0', left - 6, y(z), { align: 'right', color: MUTED, size }); }
    text(ctx, 'm', left - 6, top - 18, { align: 'right', color: MUTED, size });
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
    } else state.at = { lat: Math.round(gm.latOf(px) * 10) / 10, z: Math.round(gm.zOf(py) * 10) / 10 };
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

function trace(ctx, points, x, y, join, backward) {
  let joined = join;
  for (const [s, t] of backward ? points.slice().reverse() : points) {
    if (!Number.isFinite(s)) continue;
    if (joined) ctx.lineTo(x(s), y(t)); else ctx.moveTo(x(s), y(t));
    joined = true;
  }
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
