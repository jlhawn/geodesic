import { Figure, choice, legend, readout, text, arrow, anomalyColor, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';

const COLS = 9, ROWS = 7, SQRT3 = Math.sqrt(3);

export function mountStagger(root) {
  const controls = root.querySelector('.controls');
  const values = new Float64Array(COLS * ROWS);
  let grid = 'A', brush = 1, geometry = null, chips = [];
  const index = (c, r) => r * COLS + c;
  const inside = (c, r) => c >= 0 && c < COLS && r >= 0 && r < ROWS;
  const around = (c, r) => ((r & 1) ? [[1, 0], [-1, 0], [0, -1], [1, -1], [0, 1], [1, 1]] : [[1, 0], [-1, 0], [-1, -1], [0, -1], [-1, 1], [0, 1]]).map(([dc, dr]) => [c + dc, r + dr]);
  const patterns = {
    stripes: (c, r) => (r & 1 ? -1 : 1),
    checkerboard: (c, r) => (r & 1 ? 1 : ((c + (r >> 1)) & 1 ? -1 : 1)),
    bump: (c, r) => { const [x, y] = place(c, r, 1), [cx, cy] = place(4, 3, 1); return Math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / 4); },
    clear: () => 0,
  };
  function place(c, r, R) { return [SQRT3 * R * (c + 0.5 * (r & 1)), 1.5 * R * r]; }
  function fill(name) { for (let r = 0; r < ROWS; r++) for (let c = 0; c < COLS; c++) values[index(c, r)] = patterns[name](c, r); }

  const fig = new Figure(root, { height: 400, draw });
  legend(root, [['gradient', 'lower to higher pressure', 'cool', 'warm'], ['arrow', 'the push the grid feels', ACCENT]]);
  choice(controls, { label: 'The numbers live', options: [['all at the center (A-grid)', 'A'], ['wind on the edges (C-grid)', 'C']], value: grid, onChange: (v) => { grid = v; show(); fig.render(); }, span: true });
  choice(controls, { label: 'Pressure pattern', options: [['bump', 'bump'], ['alternating rows', 'stripes'], ['checkerboard', 'checkerboard'], ['clear', 'clear']], value: 'stripes', onChange: (v) => { fill(v); show(); fig.render(); }, span: true });
  const swatches = document.createElement('div');
  swatches.className = 'choice';
  swatches.setAttribute('role', 'group');
  swatches.style.gridColumn = '1 / -1';
  const caption = document.createElement('span');
  caption.className = 'caption';
  caption.textContent = 'Brush pressure';
  caption.id = 'brush-caption';
  swatches.setAttribute('aria-labelledby', caption.id);
  swatches.append(caption);
  chips = [-1, -0.5, 0, 0.5, 1].map((value) => {
    const chip = document.createElement('button');
    chip.type = 'button';
    chip.className = 'swatch-chip';
    chip.dataset.value = value;
    chip.title = `${value > 0 ? '+' : ''}${value}`;
    chip.setAttribute('aria-label', `pressure ${chip.title}`);
    chip.addEventListener('click', () => { brush = value; paintChips(); });
    swatches.append(chip);
    return chip;
  });
  controls.append(swatches);
  function paintChips() { for (const chip of chips) { const value = Number(chip.dataset.value); chip.style.background = anomalyColor(value, 0.75); chip.setAttribute('aria-pressed', String(value === brush)); } }
  paintChips();
  const out = readout(controls);

  function pushA(c, r) {
    let gx = 0, gy = 0;
    const [x0, y0] = place(c, r, 1);
    for (const [cc, rr] of around(c, r)) {
      if (!inside(cc, rr)) return null;
      const [x1, y1] = place(cc, rr, 1), d = Math.hypot(x1 - x0, y1 - y0), p = 0.5 * (values[index(c, r)] + values[index(cc, rr)]);
      gx += p * (x1 - x0) / d; gy += p * (y1 - y0) / d;
    }
    return [-gx / 3, -gy / 3];
  }

  function edges() {
    const list = [];
    for (let r = 0; r < ROWS; r++) for (let c = 0; c < COLS; c++) for (const [cc, rr] of around(c, r)) {
      if (!inside(cc, rr) || index(cc, rr) < index(c, r)) continue;
      list.push({ a: [c, r], b: [cc, rr], push: values[index(c, r)] - values[index(cc, rr)] });
    }
    return list;
  }

  function strongest() {
    let worst = 0;
    if (grid === 'A') { for (let r = 1; r < ROWS - 1; r++) for (let c = 1; c < COLS - 1; c++) { const g = pushA(c, r); if (g) worst = Math.max(worst, Math.hypot(...g)); } }
    else for (const e of edges()) worst = Math.max(worst, Math.abs(e.push));
    return worst;
  }

  function show() {
    const worst = strongest(), flat = values.every((v) => Math.abs(v) < 1e-9);
    out.set([['strongest push anywhere', flat ? 'none: the pressure is the same everywhere' : worst < 1e-9 ? 'none: the grid cannot feel this pattern' : worst.toFixed(2)]]);
  }

  function draw(ctx, w, h) {
    paintChips();
    const R = Math.min(30, (w - 40) / (SQRT3 * (COLS + 0.5)), (h - 40) / (1.5 * (ROWS - 1) + 2));
    const x0 = (w - SQRT3 * R * (COLS + 0.5)) / 2 + SQRT3 * R / 2, y0 = (h - 1.5 * R * (ROWS - 1)) / 2;
    const center = (c, r) => { const [x, y] = place(c, r, R); return [x0 + x, y0 + y]; };
    geometry = { R, center };
    for (let r = 0; r < ROWS; r++) for (let c = 0; c < COLS; c++) {
      const [x, y] = center(c, r);
      ctx.beginPath();
      for (let k = 0; k < 6; k++) { const a = Math.PI / 6 + k * Math.PI / 3; ctx.lineTo(x + R * Math.cos(a), y + R * Math.sin(a)); }
      ctx.closePath();
      ctx.fillStyle = anomalyColor(values[index(c, r)], 0.75); ctx.fill();
      ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.stroke();
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(x, y, 2.5, 0, Math.PI * 2); ctx.fill();
    }
    const scale = R * 1.6;
    if (grid === 'A') {
      for (let r = 0; r < ROWS; r++) for (let c = 0; c < COLS; c++) {
        const g = pushA(c, r);
        if (!g) continue;
        const [x, y] = center(c, r), len = Math.hypot(...g) * scale;
        if (len > 1.5) arrow(ctx, x - g[0] * scale / 2, y - g[1] * scale / 2, x + g[0] * scale / 2, y + g[1] * scale / 2, { color: ACCENT, width: 2, head: 7 });
      }
    } else {
      for (const e of edges()) {
        const [xa, ya] = center(...e.a), [xb, yb] = center(...e.b), mx = (xa + xb) / 2, my = (ya + yb) / 2, d = Math.hypot(xb - xa, yb - ya);
        const nx = (xb - xa) / d, ny = (yb - ya) / d;
        ctx.strokeStyle = 'rgba(255,255,255,0.8)'; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(mx - nx * 5, my - ny * 5); ctx.lineTo(mx + nx * 5, my + ny * 5); ctx.stroke();
        const len = e.push * scale * 0.5;
        if (Math.abs(len) > 1.5) arrow(ctx, mx - nx * len / 2, my - ny * len / 2, mx + nx * len / 2, my + ny * len / 2, { color: ACCENT, width: 2, head: 7 });
      }
    }
    text(ctx, grid === 'A' ? 'pressure and wind both at the dots' : 'pressure at the dots, wind on the white ticks', w / 2, h - 12, { align: 'center', color: MUTED, size: 11 });
    text(ctx, 'click or drag across cells to paint them with the brush pressure', w / 2, 12, { align: 'center', color: MUTED, size: 11 });
    void INK;
  }

  function cellAt({ x: px, y: py }) {
    if (!geometry) return null;
    for (let r = 0; r < ROWS; r++) for (let c = 0; c < COLS; c++) {
      const [x, y] = geometry.center(c, r);
      if (Math.hypot(px - x, py - y) < geometry.R * 0.9) return index(c, r);
    }
    return null;
  }
  function paint(p) {
    const i = cellAt(p);
    if (i === null || values[i] === brush) return;
    values[i] = brush;
    show();
    fig.render();
  }
  fig.pointer({ hit: (p) => cellAt(p) !== null, down: paint, move: paint });

  fill('stripes');
  show();
}
