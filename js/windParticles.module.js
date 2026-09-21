/*
 * Wind shown by particles: each frame every particle moves with the
 * field at its own position and is drawn as a dot whose opacity rises
 * with the wind speed, from `minimumOpacity` at rest to full at the
 * reference speed. Particles live a few seconds and respawn at random
 * so the field stays evenly seeded where the flow converges.
 */
export function createWindParticles(container, viewer, grid, { density = 0.02, referenceSpeed = 15, pixelsPerFrame = 0.5, size = 2, minimumOpacity = 0.25 } = {}) {
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
  let width = 0, height = 0, dpr = 1, count = 0;
  const capacity = 400000;
  const position = new Float32Array(3 * capacity);
  const cellOf = new Int32Array(capacity);
  const age = new Uint16Array(capacity);
  const lifetime = new Uint16Array(capacity);
  let field = null, reference = referenceSpeed;

  function resize() {
    dpr = window.devicePixelRatio || 1;
    width = container.clientWidth; height = container.clientHeight;
    canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
    context.setTransform(dpr, 0, 0, dpr, 0, 0);
    const wanted = Math.min(capacity, Math.round(density * width * height));
    for (let n = count; n < wanted; n++) { spawn(n); age[n] = Math.floor(Math.random() * lifetime[n]); }
    count = wanted;
  }

  function spawn(n) {
    const i = Math.floor(Math.random() * C);
    const cx = centers[3 * i], cy = centers[3 * i + 1], cz = centers[3 * i + 2];
    let jx = Math.random() - 0.5, jy = Math.random() - 0.5, jz = Math.random() - 0.5;
    const radial = jx * cx + jy * cy + jz * cz;
    jx -= radial * cx; jy -= radial * cy; jz -= radial * cz;
    let x = cx + cellSpan * jx, y = cy + cellSpan * jy, z = cz + cellSpan * jz;
    const norm = Math.hypot(x, y, z);
    position[3 * n] = x / norm; position[3 * n + 1] = y / norm; position[3 * n + 2] = z / norm;
    cellOf[n] = i;
    age[n] = 0;
    lifetime[n] = 60 + Math.floor(Math.random() * 180);
  }
  function reset() {
    context.clearRect(0, 0, width, height);
    for (let n = 0; n < count; n++) { spawn(n); age[n] = Math.floor(Math.random() * lifetime[n]); }
  }
  const resizeObserver = new ResizeObserver(resize);
  resizeObserver.observe(container);
  resize();

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

  function frame() {
    if (!running) return;
    requestAnimationFrame(frame);
    if (!field || width === 0) return;
    const step = pixelsPerFrame / (reference * viewer.pixelsPerUnit());
    const paths = Array.from({ length: buckets }, () => new Path2D());
    for (let n = 0; n < count; n++) {
      if (++age[n] >= lifetime[n]) { spawn(n); continue; }
      const v = sample(n);
      const speed = Math.hypot(v[0], v[1], v[2]);
      let x = position[3 * n] + step * v[0], y = position[3 * n + 1] + step * v[1], z = position[3 * n + 2] + step * v[2];
      const norm = Math.hypot(x, y, z);
      x /= norm; y /= norm; z /= norm;
      position[3 * n] = x; position[3 * n + 1] = y; position[3 * n + 2] = z;
      viewer.projectPoint(x, y, z, screen);
      if (screen[2] <= 0) continue;
      const bucket = Math.round((buckets - 1) * Math.min(1, speed / reference));
      paths[bucket].rect(screen[0] - size / 2, screen[1] - size / 2, size, size);
    }
    context.clearRect(0, 0, width, height);
    for (let b = 0; b < buckets; b++) {
      context.fillStyle = `rgba(255, 255, 255, ${opacity[b]})`;
      context.fill(paths[b]);
    }
  }
  requestAnimationFrame(frame);

  return {
    setField(vectors, speed = referenceSpeed) { field = vectors; reference = speed; },
    setVisible(visible) { canvas.style.display = visible ? 'block' : 'none'; if (!visible) context.clearRect(0, 0, width, height); },
    reset,
    dispose() { running = false; resizeObserver.disconnect(); canvas.remove(); },
  };
}
