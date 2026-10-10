import { R, CP, KAPPA, G, P0 } from './physics.module.js';
export { R, CP, KAPPA, G, P0 };

export const ACCENT = '#ffe8a0', INK = '#ddd', MUTED = '#8a8a8a', GRID = 'rgba(255,255,255,0.1)', LINE = 'rgba(255,255,255,0.35)';
export const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
export const lerp = (a, b, t) => a + (b - a) * t;

const figures = new Map();
const observer = new IntersectionObserver((entries) => { for (const entry of entries) figures.get(entry.target)?.setVisible(entry.isIntersecting); }, { rootMargin: '120px' });

export class Figure {
  constructor(root, { height = 360, minHeight = 260, step = null, draw, context = '2d' }) {
    this.root = root;
    this.stage = root.querySelector('.stage');
    this.canvas = document.createElement('canvas');
    this.stage.append(this.canvas);
    this.ctx = context === '2d' ? this.canvas.getContext('2d') : null;
    Object.assign(this, { height, minHeight, step, draw, visible: false, running: false, frame: 0, last: 0, width: 0, h: 0, dpr: 1 });
    this.tick = this.tick.bind(this);
    new ResizeObserver(() => this.resize()).observe(this.stage);
    figures.set(root, this);
    observer.observe(root);
    this.resize();
  }

  resize() {
    const w = this.stage.clientWidth;
    if (!w) return;
    const h = Math.round(Math.max(this.minHeight, Math.min(this.height, w * 0.8)));
    this.dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(w * this.dpr);
    this.canvas.height = Math.round(h * this.dpr);
    this.canvas.style.height = `${h}px`;
    this.width = w;
    this.h = h;
    matchMedia(`(resolution: ${this.dpr}dppx)`).addEventListener('change', () => this.resize(), { once: true });
    this.render();
  }

  setVisible(on) { this.visible = on; this.sync(); }
  play(on = true) { this.running = on; this.sync(); }

  sync() {
    const should = this.visible && this.running && this.step;
    if (should && !this.frame) { this.last = performance.now(); this.frame = requestAnimationFrame(this.tick); }
    else if (!should && this.frame) { cancelAnimationFrame(this.frame); this.frame = 0; }
  }

  tick(now) {
    this.frame = 0;
    const dt = clamp((now - this.last) / 1000, 0, 0.25);
    this.last = now;
    this.step(dt);
    this.render();
    this.sync();
  }

  render() {
    if (!this.width) return;
    const { ctx, width: w, h, dpr } = this;
    if (ctx) { ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h); }
    this.draw(ctx, w, h);
  }

  pointer({ hit = () => true, down, move, up }) {
    const at = (e) => { const r = this.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
    this.canvas.addEventListener('touchstart', (e) => { if (hit(at(e.touches[0]))) e.preventDefault(); }, { passive: false });
    this.canvas.addEventListener('pointerdown', (e) => { if (hit(at(e)) && down?.(at(e), e) !== false) { this.canvas.setPointerCapture(e.pointerId); e.preventDefault(); } });
    this.canvas.addEventListener('pointermove', (e) => { if (this.canvas.hasPointerCapture(e.pointerId)) move?.(at(e), e); });
    const release = (e) => { if (this.canvas.hasPointerCapture(e.pointerId)) { this.canvas.releasePointerCapture(e.pointerId); up?.(at(e), e); } };
    this.canvas.addEventListener('pointerup', release);
    this.canvas.addEventListener('pointercancel', release);
  }
}

const el = (tag, cls, text) => { const e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; };
let ids = 0;

export function slider(parent, { label, min, max, step = 'any', value, format = (v) => `${v}`, onInput = () => {}, span = false }) {
  const wrap = el('div', 'control'), lab = el('label', null, label), input = el('input'), out = el('output');
  if (span) wrap.style.gridColumn = '1 / -1';
  Object.assign(input, { type: 'range', min, max, step, value, id: `control${++ids}` });
  lab.htmlFor = input.id;
  const show = () => { out.textContent = format(Number(input.value)); };
  input.addEventListener('input', () => { show(); onInput(Number(input.value)); });
  wrap.append(lab, input, out);
  parent.append(wrap);
  show();
  return { get value() { return Number(input.value); }, set value(v) { input.value = v; show(); } };
}

export function choice(parent, { label, options, value, onChange = () => {}, span = false }) {
  const wrap = el('div', 'choice');
  if (span) wrap.style.gridColumn = '1 / -1';
  wrap.setAttribute('role', 'group');
  if (label) { const caption = el('span', 'caption', label); caption.id = `caption${++ids}`; wrap.setAttribute('aria-labelledby', caption.id); wrap.append(caption); }
  let current = value;
  const set = (v) => { current = v; for (const b of made) b.setAttribute('aria-pressed', String(b.dataset.value === v)); };
  const made = options.map(([text, v]) => { const b = el('button', null, text); b.type = 'button'; b.dataset.value = v; b.addEventListener('click', () => { set(v); onChange(v); }); wrap.append(b); return b; });
  set(value);
  parent.append(wrap);
  return { get value() { return current; }, set value(v) { set(v); } };
}

export function buttons(parent, items) {
  const wrap = el('div', 'buttons');
  const made = items.map(([text, fn]) => { const b = el('button', null, text); b.type = 'button'; b.addEventListener('click', () => fn(b)); wrap.append(b); return b; });
  parent.append(wrap);
  return made;
}

export function readout(parent) {
  const box = el('div', 'readout');
  parent.append(box);
  let frames = 0;
  return {
    set(pairs, every = 1) {
      if (frames++ % every) return;
      box.replaceChildren(...pairs.flatMap(([key, value], i) => [i ? ` · ${key} ` : `${key} `, el('b', null, value)]));
    },
  };
}

const SWATCHES = {
  gradient: (a, b) => `<span class="gradient" style="background: linear-gradient(90deg, ${a}, transparent 50%, ${b})"></span>`,
  ramp: (a, b, mid = null) => `<span class="gradient" style="background: linear-gradient(90deg, ${a}, ${mid ?? b}, ${b})"></span>`,
  force: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 6h19" stroke="${color}" stroke-width="1.5" stroke-dasharray="3 2.5" stroke-linecap="round"/><path d="M28 6l-6-3v6z" fill="none" stroke="${color}" stroke-width="1.5" stroke-linejoin="round"/></svg>`,
  arrow: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 6h22" stroke="${color}" stroke-width="1.5" stroke-linecap="round"/><path d="M28 6l-6-3v6z" fill="${color}"/></svg>`,
  varrow: (color) => `<svg width="30" height="16" viewBox="0 0 30 16"><path d="M10 14V4" stroke="${color}" stroke-width="1.5" stroke-linecap="round"/><path d="M10 1l-3 6h6z" fill="${color}"/><path d="M20 2v10" stroke="${color}" stroke-width="1.5" stroke-linecap="round"/><path d="M20 15l-3-6h6z" fill="${color}"/></svg>`,
  bar: (color) => `<svg width="30" height="14" viewBox="0 0 30 14"><path d="M2 4h26" stroke="rgba(255,255,255,0.35)"/><rect x="9" y="4" width="12" height="8" fill="${color}"/></svg>`,
  line: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 9c6 0 6-6 12-6s6 6 12 6" fill="none" stroke="${color}" stroke-width="2"/></svg>`,
  dash: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 6h26" stroke="${color}" stroke-width="1.5" stroke-dasharray="3 4"/></svg>`,
  faint: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 6h26" stroke="${color}" stroke-width="1"/></svg>`,
  dots: (color) => `<svg width="30" height="12" viewBox="0 0 30 12"><path d="M2 6h26" stroke="${color}" stroke-width="2" stroke-dasharray="1.5 3.5" stroke-linecap="round"/></svg>`,
};

const legends = [];
const resolve = (arg) => typeof arg === 'string' && arg in palette ? rgba(palette[arg], 0.85) : arg;

function fillLegend(wrap, items) {
  wrap.replaceChildren();
  for (const [kind, label, ...args] of items) {
    const item = el('span', 'item'), swatch = el('span', 'swatch');
    swatch.innerHTML = SWATCHES[kind](...args.map(resolve));
    item.append(swatch, document.createTextNode(label));
    wrap.append(item);
  }
}

export function caption(root, content) {
  const line = el('div', 'caption', content);
  root.querySelector('.stage').before(line);
  return line;
}

export function legend(root, items) {
  const wrap = el('div', 'legend'), entry = [wrap, items];
  root.querySelector('.stage').after(wrap);
  legends.push(entry);
  fillLegend(wrap, items);
  return { set(next) { entry[1] = next; fillLegend(wrap, next); } };
}

export function arrow(ctx, x0, y0, x1, y1, { color = INK, width = 1.5, head = 6, dash = null, open = false } = {}) {
  const dx = x1 - x0, dy = y1 - y0, len = Math.hypot(dx, dy);
  if (len < 0.5) return;
  const ux = dx / len, uy = dy / len, hd = Math.min(head, len);
  ctx.save();
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = width; ctx.lineCap = 'round'; ctx.lineJoin = 'round';
  if (dash) ctx.setLineDash(dash);
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1 - ux * hd * (open ? 1 : 0.7), y1 - uy * hd * (open ? 1 : 0.7)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - ux * hd - uy * hd * 0.5, y1 - uy * hd + ux * hd * 0.5); ctx.lineTo(x1 - ux * hd + uy * hd * 0.5, y1 - uy * hd - ux * hd * 0.5); ctx.closePath();
  if (open) ctx.stroke(); else ctx.fill();
  ctx.restore();
}

export function text(ctx, str, x, y, { color = INK, size = 12, align = 'left', baseline = 'middle', weight = 400 } = {}) {
  ctx.save();
  ctx.fillStyle = color; ctx.font = `${weight} ${size}px system-ui, sans-serif`; ctx.textAlign = align; ctx.textBaseline = baseline;
  ctx.fillText(str, x, y);
  ctx.restore();
}

const PALETTES = {
  standard: { cool: [80, 150, 255], warm: [255, 110, 70], neutral: [150, 150, 150] },
  safe: { cool: [0, 114, 178], warm: [230, 159, 0], neutral: [150, 150, 150] },
};
let palette = PALETTES.standard;
export let paletteVersion = 0;
const rgba = (c, a = 1) => `rgba(${c[0]}, ${c[1]}, ${c[2]}, ${a})`;

export function anomalyColor(t, alpha = 1) {
  return rgba(t >= 0 ? palette.warm : palette.cool, clamp(Math.abs(t), 0, 1) * alpha);
}

export function rampColor(f) {
  const stops = [palette.cool, palette.neutral, palette.warm];
  const t = clamp(f, 0, 1) * 2, i = Math.min(1, Math.floor(t)), u = t - i;
  const [a, b] = [stops[i], stops[i + 1]];
  return `rgb(${Math.round(lerp(a[0], b[0], u))}, ${Math.round(lerp(a[1], b[1], u))}, ${Math.round(lerp(a[2], b[2], u))})`;
}

export const DARK_NEUTRAL = [52, 55, 62];

export function rampRGB(f, out = [0, 0, 0], neutral = palette.neutral) {
  const stops = [palette.cool, neutral, palette.warm];
  const t = clamp(f, 0, 1) * 2, i = Math.min(1, Math.floor(t)), u = t - i, a = stops[i], b = stops[i + 1];
  out[0] = lerp(a[0], b[0], u) / 255; out[1] = lerp(a[1], b[1], u) / 255; out[2] = lerp(a[2], b[2], u) / 255;
  return out;
}

export function setPalette(name) {
  palette = PALETTES[name] ?? PALETTES.standard;
  paletteVersion++;
  for (const [wrap, items] of legends) fillLegend(wrap, items);
  for (const figure of figures.values()) figure.render();
  try { localStorage.setItem('explain-palette', name); } catch {}
}

export function paletteControl(parent) {
  let saved = 'standard';
  try { saved = localStorage.getItem('explain-palette') ?? 'standard'; } catch {}
  if (!(saved in PALETTES)) saved = 'standard';
  if (saved !== 'standard') setPalette(saved);
  choice(parent, { label: 'Colors', options: [['standard', 'standard'], ['color-blind safe', 'safe']], value: saved, onChange: setPalette });
}

export function thermometer(ctx, x, top, height, value, { min, max, label, unit = '°C', color = ACCENT, width = 10 }) {
  const bulb = width * 0.9, tubeBottom = top + height, f = clamp((value - min) / (max - min), 0, 1);
  ctx.save();
  ctx.lineWidth = 1.5; ctx.strokeStyle = LINE; ctx.fillStyle = 'rgba(255,255,255,0.06)';
  ctx.beginPath(); ctx.roundRect(x - width / 2, top, width, height, width / 2); ctx.fill(); ctx.stroke();
  ctx.fillStyle = color;
  ctx.beginPath(); ctx.roundRect(x - width / 4, tubeBottom - f * (height - width / 2) - width / 2, width / 2, f * (height - width / 2) + width / 2, width / 4); ctx.fill();
  ctx.beginPath(); ctx.arc(x, tubeBottom + bulb * 0.6, bulb, 0, Math.PI * 2); ctx.fill();
  for (let v = min + 20; v <= max; v += 20) { const y = tubeBottom - width / 4 - ((v - min) / (max - min)) * (height - width / 2); ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x + width / 2 + 2, y); ctx.lineTo(x + width / 2 + 6, y); ctx.stroke(); text(ctx, `${v}`, x + width / 2 + 9, y, { color: MUTED, size: 10 }); }
  ctx.restore();
  text(ctx, `${value.toFixed(1)} ${unit}`, x, top - 12, { align: 'center', color: '#fff', weight: 500 });
  (Array.isArray(label) ? label : [label]).forEach((line, n) => { if (line) text(ctx, line, x, tubeBottom + bulb * 2.4 + n * 12, { align: 'center', color: MUTED, size: 11 }); });
}
