/*
 * Wind shown by particles: each frame every particle moves with the
 * field at its own position and is drawn as a dot whose opacity rises
 * with the wind speed, from `minimumOpacity` at rest to full at the
 * reference speed. The screen is kept evenly covered: it is divided
 * into bins, a bin with too few particles receives new ones at random
 * points inside it, and a bin with too many loses one, so the flow
 * neither piles particles up where it converges nor empties them where
 * it diverges. Particles out of view for a second are retired,
 * and with an admission mask (the sea, for currents) particles exist
 * only on admitted cells.
 */
export function createWindParticles(container, viewer, grid, { density = 0.02, referenceSpeed = 15, pixelsPerFrame = 0.5, size = 2, minimumOpacity = 0.25, bin = 32, slack = 0.4, maximum = 200000 } = {}) {
  const C = grid.size;
  const centers = new Float32Array(3 * C);
  const neighborCount = new Uint8Array(C);
  const neighbors = new Int32Array(6 * C);
  for (const cell of grid) {
    const i = cell.index, c = cell.centerVertex;
    centers[3 * i] = c.x; centers[3 * i + 1] = c.y; centers[3 * i + 2] = c.z;
    neighborCount[i] = cell.neighbors.length;
    for (let k = 0; k < cell.neighbors.length; k++) neighbors[6 * i + k] = cell.neighbors[k].index;
  }
  const cellSpan = Math.sqrt(4 * Math.PI / C);

  const canvas = document.createElement('canvas');
  Object.assign(canvas.style, { position: 'absolute', inset: '0', width: '100%', height: '100%', pointerEvents: 'none', zIndex: '5' });
  container.appendChild(canvas);
  const context = canvas.getContext('2d');
  let width = 0, height = 0, dpr = 1;
  const position = new Float32Array(3 * maximum);
  const cellOf = new Int32Array(maximum);
  const hidden = new Uint16Array(maximum);
  const active = new Int32Array(maximum), slotOf = new Int32Array(maximum), free = new Int32Array(maximum);
  let count = 0, freeCount = 0;
  let cols = 0, rows = 0, counts = null, seed = null, seen = null, coverage = null;
  let field = null, reference = referenceSpeed, mask = null, visible = true, lastCell = 0;
  let frames = 0, coverageVersion = -1, coverageFrame = -1;

  function resize() {
    dpr = window.devicePixelRatio || 1;
    width = container.clientWidth; height = container.clientHeight;
    canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
    context.setTransform(dpr, 0, 0, dpr, 0, 0);
    cols = Math.ceil(width / bin); rows = Math.ceil(height / bin);
    counts = new Int32Array(cols * rows); seen = new Int32Array(cols * rows); seed = new Int32Array(cols * rows).fill(-1); coverage = new Float32Array(cols * rows);
    coverageVersion = -1;
  }
  function reset() {
    context.clearRect(0, 0, width, height);
    count = 0; freeCount = maximum;
    for (let i = 0; i < maximum; i++) free[i] = maximum - 1 - i;
    if (seed) seed.fill(-1);
    coverageVersion = -1;
  }
  reset();
  const resizeObserver = new ResizeObserver(resize);
  resizeObserver.observe(container);
  resize();

  const point = [0, 0, 0];
  function measureCoverage() {
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        let hits = 0;
        for (let k = 0; k < 4; k++) {
          const px = Math.min(width - 0.5, (c + (k & 1 ? 0.9 : 0.1)) * bin), py = Math.min(height - 0.5, (r + (k & 2 ? 0.9 : 0.1)) * bin);
          if (!viewer.unprojectPoint(px, py, point)) continue;
          if (!mask) { hits++; continue; }
          lastCell = locate(point[0], point[1], point[2], lastCell);
          if (mask[lastCell]) hits++;
        }
        coverage[r * cols + c] = hits / 4;
      }
    }
  }

  const dotWith = (x, y, z, j) => x * centers[3 * j] + y * centers[3 * j + 1] + z * centers[3 * j + 2];
  function locate(x, y, z, start) {
    let here = start, best = dotWith(x, y, z, here);
    for (;;) {
      let next = here;
      for (let k = 0; k < neighborCount[here]; k++) { const j = neighbors[6 * here + k]; const d = dotWith(x, y, z, j); if (d > best) { best = d; next = j; } }
      if (next === here) return here;
      here = next;
    }
  }
  function add(x, y, z, cell) {
    if (freeCount === 0) return;
    const n = free[--freeCount];
    position[3 * n] = x; position[3 * n + 1] = y; position[3 * n + 2] = z;
    cellOf[n] = cell; hidden[n] = 0;
    slotOf[n] = count; active[count++] = n;
  }
  function remove(n) {
    const a = slotOf[n], last = active[--count];
    active[a] = last; slotOf[last] = a;
    free[freeCount++] = n;
  }

  const wind = [0, 0, 0];
  function sample(n) {
    const x = position[3 * n], y = position[3 * n + 1], z = position[3 * n + 2];
    let i = cellOf[n], best = i, bestDot = -2;
    const candidates = neighborCount[i] + 1;
    let wx = 0, wy = 0, wz = 0, total = 0;
    for (let k = 0; k < candidates; k++) {
      const j = k === 0 ? i : neighbors[6 * i + k - 1];
      const dot = x * centers[3 * j] + y * centers[3 * j + 1] + z * centers[3 * j + 2];
      if (dot > bestDot) { bestDot = dot; best = j; }
      const w = 1 / (2 * (1 - dot) + 1e-4 * cellSpan * cellSpan);
      wx += w * field[3 * j]; wy += w * field[3 * j + 1]; wz += w * field[3 * j + 2];
      total += w;
    }
    cellOf[n] = best;
    wind[0] = wx / total; wind[1] = wy / total; wind[2] = wz / total;
    return wind;
  }

  const screen = [0, 0, 0];
  const buckets = 8;
  const opacity = Array.from({ length: buckets }, (_, b) => (minimumOpacity + (1 - minimumOpacity) * b / (buckets - 1)).toFixed(3));
  let running = true;
  const retiring = [];

  function frame() {
    if (!running) return;
    requestAnimationFrame(frame);
    if (!field || !visible || width === 0) return;
    frames++;
    const version = viewer.viewVersion();
    if (version !== coverageVersion && (coverageVersion === -1 || frames - coverageFrame >= 6)) { measureCoverage(); coverageVersion = version; coverageFrame = frames; }
    const step = pixelsPerFrame / (reference * viewer.pixelsPerUnit());
    const paths = Array.from({ length: buckets }, () => new Path2D());
    counts.fill(0); seen.fill(-1);
    retiring.length = 0;
    for (let a = 0; a < count; a++) {
      const n = active[a];
      const v = sample(n);
      if (mask && !mask[cellOf[n]]) { retiring.push(n); continue; }
      const speed = Math.hypot(v[0], v[1], v[2]);
      let x = position[3 * n] + step * v[0], y = position[3 * n + 1] + step * v[1], z = position[3 * n + 2] + step * v[2];
      const norm = Math.hypot(x, y, z);
      x /= norm; y /= norm; z /= norm;
      position[3 * n] = x; position[3 * n + 1] = y; position[3 * n + 2] = z;
      viewer.projectPoint(x, y, z, screen);
      if (screen[2] <= 0 || screen[0] < 0 || screen[0] >= width || screen[1] < 0 || screen[1] >= height) { if (++hidden[n] > 60) retiring.push(n); continue; }
      hidden[n] = 0;
      const b = (screen[1] / bin | 0) * cols + (screen[0] / bin | 0);
      counts[b]++; seen[b] = n; seed[b] = cellOf[n];
      const bucket = Math.round((buckets - 1) * Math.min(1, speed / reference));
      paths[bucket].rect(screen[0] - size / 2, screen[1] - size / 2, size, size);
    }
    for (const n of retiring) remove(n);

    let budget = 2000;
    const area = bin * bin;
    for (let b = 0; b < counts.length && budget > 0; b++) {
      const target = density * area * coverage[b];
      if (target < 0.5) continue;
      if (counts[b] > target * (1 + slack)) { if (seen[b] >= 0) remove(seen[b]); continue; }
      if (counts[b] >= target * (1 - slack)) continue;
      const c = b % cols, r = (b - c) / cols;
      for (let wanted = Math.min(4, Math.ceil(target - counts[b])); wanted > 0 && budget > 0; wanted--) {
        const px = (c + Math.random()) * bin, py = (r + Math.random()) * bin;
        if (px >= width || py >= height || !viewer.unprojectPoint(px, py, point)) continue;
        const cell = locate(point[0], point[1], point[2], seed[b] >= 0 ? seed[b] : lastCell);
        lastCell = cell; seed[b] = cell;
        if (mask && !mask[cell]) continue;
        add(point[0], point[1], point[2], cell);
        budget--;
      }
    }

    context.clearRect(0, 0, width, height);
    for (let b = 0; b < buckets; b++) {
      context.fillStyle = `rgba(255, 255, 255, ${opacity[b]})`;
      context.fill(paths[b]);
    }
  }
  requestAnimationFrame(frame);

  return {
    setField(vectors, speed = referenceSpeed, admit = null) { field = vectors; reference = speed; if (admit !== mask) { mask = admit; coverageVersion = -1; } },
    setVisible(show) { visible = show; canvas.style.display = show ? 'block' : 'none'; if (!show) context.clearRect(0, 0, width, height); },
    reset,
    count: () => count,
    dispose() { running = false; resizeObserver.disconnect(); canvas.remove(); },
  };
}
