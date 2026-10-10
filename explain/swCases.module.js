export const EARTH = { a: 6.37122e6, omega: 7.292e-5, g: 9.80616 };

export const lonLat = (x, y, z) => ({ lon: Math.atan2(y, x), lat: Math.atan2(z, Math.hypot(x, y)) });

export function edgeNormalVelocity(mesh, windAt) {
  const { nEdges, xEdge, nEdge } = mesh, u = new Float64Array(nEdges);
  for (let e = 0; e < nEdges; e++) {
    const { lon, lat } = lonLat(xEdge[3 * e], xEdge[3 * e + 1], xEdge[3 * e + 2]);
    const { zonal, meridional } = windAt(lon, lat);
    const cosLon = Math.cos(lon), sinLon = Math.sin(lon), cosLat = Math.cos(lat), sinLat = Math.sin(lat);
    const vx = -sinLon * zonal - sinLat * cosLon * meridional, vy = cosLon * zonal - sinLat * sinLon * meridional, vz = cosLat * meridional;
    u[e] = vx * nEdge[3 * e] + vy * nEdge[3 * e + 1] + vz * nEdge[3 * e + 2];
  }
  return u;
}

export function cellField(mesh, valueAt) {
  const { nCells, xCell } = mesh, h = new Float64Array(nCells);
  for (let i = 0; i < nCells; i++) { const { lon, lat } = lonLat(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]); h[i] = valueAt(lon, lat); }
  return h;
}

export function hyperdiffusion(mesh, hours) {
  let sum = 0;
  for (let e = 0; e < mesh.nEdges; e++) sum += mesh.dcEdge[e];
  return (sum / mesh.nEdges / Math.PI) ** 4 / (hours * 3600);
}

export function bump({ depth = 4000, height = 100, width = 1.5e6, lat0 = Math.PI / 4, lon0 = 0 } = {}) {
  const c0 = [Math.cos(lat0) * Math.cos(lon0), Math.cos(lat0) * Math.sin(lon0), Math.sin(lat0)];
  return {
    mean: depth,
    height: (lon, lat) => {
      const p = [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
      const d = EARTH.a * Math.acos(Math.min(1, p[0] * c0[0] + p[1] * c0[1] + p[2] * c0[2]));
      return depth + height * Math.exp(-((d / width) ** 2));
    },
    wind: () => ({ zonal: 0, meridional: 0 }),
  };
}

export function galewsky({ perturbed = true } = {}) {
  const { a, omega, g } = EARTH;
  const phi0 = Math.PI / 7, phi1 = Math.PI / 2 - phi0, uMax = 80, en = Math.exp(-4 / (phi1 - phi0) ** 2), hMean = 1e4;
  const jet = (lat) => (lat <= phi0 || lat >= phi1 ? 0 : (uMax / en) * Math.exp(1 / ((lat - phi0) * (lat - phi1))));
  const steps = 20000, dphi = Math.PI / steps, table = new Float64Array(steps + 1);
  let integral = 0;
  for (let k = 1; k <= steps; k++) {
    const mid = -Math.PI / 2 + (k - 0.5) * dphi, u = jet(mid);
    integral += a * u * (2 * omega * Math.sin(mid) + Math.tan(mid) * u / a) * dphi;
    table[k] = -integral / g;
  }
  let weighted = 0, area = 0;
  for (let k = 0; k <= steps; k++) { const w = Math.cos(-Math.PI / 2 + k * dphi); weighted += w * table[k]; area += w; }
  const offset = hMean - weighted / area;
  const balanced = (lat) => { const x = (lat + Math.PI / 2) / dphi, k = Math.min(steps - 1, Math.max(0, Math.floor(x))), f = x - k; return offset + table[k] * (1 - f) + table[k + 1] * f; };
  const bumpAt = (lon, lat) => 120 * Math.cos(lat) * Math.exp(-((lon / (1 / 3)) ** 2)) * Math.exp(-(((Math.PI / 4 - lat) / (1 / 15)) ** 2));
  return { mean: hMean, height: (lon, lat) => balanced(lat) + (perturbed ? bumpAt(lon, lat) : 0), wind: (lon, lat) => ({ zonal: jet(lat), meridional: 0 }) };
}

export function cellVelocity(mesh, u, i, out, k, offset = 0) {
  const { maxEdges, nEdgesOnCell, edgesOnCell, nEdge, xCell } = mesh;
  const x = xCell[3 * i], y = xCell[3 * i + 1], z = xCell[3 * i + 2];
  const rho = Math.hypot(x, y) || 1e-12, east = [-y / rho, x / rho, 0], north = [-z * x / rho, -z * y / rho, rho];
  let saa = 0, sab = 0, sbb = 0, sau = 0, sbu = 0;
  for (let m = 0; m < nEdgesOnCell[i]; m++) {
    const e = edgesOnCell[maxEdges * i + m], n0 = nEdge[3 * e], n1 = nEdge[3 * e + 1], n2 = nEdge[3 * e + 2];
    const a = n0 * east[0] + n1 * east[1] + n2 * east[2], b = n0 * north[0] + n1 * north[1] + n2 * north[2], ue = u[offset + e];
    saa += a * a; sab += a * b; sbb += b * b; sau += a * ue; sbu += b * ue;
  }
  const det = saa * sbb - sab * sab, ue = (sau * sbb - sbu * sab) / det, un = (saa * sbu - sab * sau) / det;
  if (out.length === 2) { out[0] = ue; out[1] = un; return out; }
  for (let c = 0; c < 3; c++) out[3 * k + c] = ue * east[c] + un * north[c];
  return out;
}

export function haurwitz() {
  const { a, omega, g } = EARTH, w = 7.848e-6, K = 7.848e-6, R = 4, h0 = 8e3;
  return {
    mean: h0,
    wind: (lon, lat) => {
      const c = Math.cos(lat), s = Math.sin(lat);
      return { zonal: a * w * c + a * K * c ** (R - 1) * (R * s * s - c * c) * Math.cos(R * lon), meridional: -a * K * R * c ** (R - 1) * s * Math.sin(R * lon) };
    },
    height: (lon, lat) => {
      const c = Math.cos(lat);
      const A = (w / 2) * (2 * omega + w) * c * c + (K * K / 4) * c ** (2 * R) * ((R + 1) * c * c + (2 * R * R - R - 2) - 2 * R * R / (c * c));
      const B = (2 * (omega + w) * K / ((R + 1) * (R + 2))) * c ** R * ((R * R + 2 * R + 2) - (R + 1) ** 2 * c * c);
      const C = (K * K / 4) * c ** (2 * R) * ((R + 1) * c * c - (R + 2));
      return h0 + (a * a / g) * (A + B * Math.cos(R * lon) + C * Math.cos(2 * R * lon));
    },
  };
}
