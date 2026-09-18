/*
 * Wind traced by particles: each particle rides the wind field over the
 * sphere and leaves a fading trail on a canvas laid over the globe, the
 * way earth.nullschool.net draws wind. Trails are brighter for faster
 * wind. The trails live in screen space, so they are cleared whenever
 * the view moves; the particles themselves stay on the globe.
 */
export function createWindParticles(container, viewer, grid, { density = 0.006, fade = 0.993, referenceSpeed = 15, pixelsPerFrame = 0.1875 } = {}) {
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
  const capacity = 100000;
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

  const from = [0, 0, 0], to = [0, 0, 0];
  const buckets = 8;
  const paths = Array.from({ length: buckets }, () => new Path2D());
  let lastVersion = -1, running = true, frames = 0;
  // Fading by a fraction of a percent per frame leaves a permanent haze
  // because 8-bit alpha rounds back to itself; fading every twelfth frame
  // by the compounded factor takes the same time to fade but reaches near
  // zero.
  const fadeEvery = 12, fadeStep = fade ** fadeEvery;

  function frame() {
    if (!running) return;
    requestAnimationFrame(frame);
    if (!field || width === 0) return;
    if (viewer.viewVersion() !== lastVersion) { lastVersion = viewer.viewVersion(); context.clearRect(0, 0, width, height); }
    if (++frames % fadeEvery === 0) {
      context.globalCompositeOperation = 'destination-in';
      context.fillStyle = `rgba(0, 0, 0, ${fadeStep})`;
      context.fillRect(0, 0, width, height);
      context.globalCompositeOperation = 'source-over';
    }

    const step = pixelsPerFrame / (reference * viewer.pixelsPerUnit());
    for (let b = 0; b < buckets; b++) paths[b] = new Path2D();
    for (let n = 0; n < count; n++) {
      if (++age[n] >= lifetime[n]) { spawn(n); continue; }
      const v = sample(n);
      const speed = Math.hypot(v[0], v[1], v[2]);
      const x = position[3 * n], y = position[3 * n + 1], z = position[3 * n + 2];
      viewer.projectPoint(x, y, z, from);
      let nx = x + step * v[0], ny = y + step * v[1], nz = z + step * v[2];
      const norm = Math.hypot(nx, ny, nz);
      nx /= norm; ny /= norm; nz /= norm;
      position[3 * n] = nx; position[3 * n + 1] = ny; position[3 * n + 2] = nz;
      viewer.projectPoint(nx, ny, nz, to);
      if (from[2] <= 0 || to[2] <= 0) continue;
      if (Math.abs(to[0] - from[0]) + Math.abs(to[1] - from[1]) > 40) continue;
      const bucket = Math.min(buckets - 1, Math.floor(buckets * Math.min(0.999, speed / reference)));
      paths[bucket].moveTo(from[0], from[1]);
      paths[bucket].lineTo(to[0], to[1]);
    }
    context.lineWidth = 1.2;
    for (let b = 0; b < buckets; b++) {
      context.strokeStyle = `rgba(255, 255, 255, ${(0.2 + 0.8 * (b + 0.5) / buckets).toFixed(3)})`;
      context.stroke(paths[b]);
    }
  }
  requestAnimationFrame(frame);

  return {
    setField(vectors, speed = referenceSpeed) { field = vectors; reference = speed; },
    setVisible(visible) { canvas.style.display = visible ? 'block' : 'none'; if (!visible) context.clearRect(0, 0, width, height); },
    dispose() { running = false; resizeObserver.disconnect(); canvas.remove(); },
  };
}
