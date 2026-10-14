import { createHash } from 'node:crypto';

const HERE = new URL('../../', import.meta.url).href;

const generic = () => new Proxy(function () {}, { get: (t, p) => (p === 'then' ? undefined : p === Symbol.toPrimitive ? () => 0 : generic()), apply: () => generic() });
const buffer = ({ size = 0 } = {}) => ({ size, label: '', mapAsync: async () => {}, getMappedRange: () => new ArrayBuffer(size), unmap() {}, destroy() {} });
const encoder = () => ({ beginComputePass: () => ({ setPipeline() {}, setBindGroup() {}, dispatchWorkgroups() {}, dispatchWorkgroupsIndirect() {}, end() {}, writeTimestamp() {} }), copyBufferToBuffer() {}, clearBuffer() {}, resolveQuerySet() {}, writeTimestamp() {}, finish: () => ({}) });
const limits = { maxStorageBuffersPerShaderStage: 10, maxStorageBufferBindingSize: 2 ** 31 - 4, maxBufferSize: 2 ** 32, maxComputeWorkgroupsPerDimension: 65535, maxComputeWorkgroupSizeX: 256, maxComputeInvocationsPerWorkgroup: 256 };

let modules = [];
const device = new Proxy({
  limits, features: new Set(), lost: new Promise(() => {}), addEventListener() {}, destroy() {},
  queue: { writeBuffer() {}, submit() {}, onSubmittedWorkDone: async () => {} },
  createBuffer: buffer, createCommandEncoder: encoder,
  createShaderModule: ({ code, label }) => { modules.push({ label, code }); return { getCompilationInfo: async () => ({ messages: [] }) }; },
  createComputePipeline: () => ({ getBindGroupLayout: () => ({}) }), createComputePipelineAsync: async () => ({ getBindGroupLayout: () => ({}) }),
  createBindGroupLayout: () => ({}), createPipelineLayout: () => ({}), createBindGroup: () => ({}), createQuerySet: () => ({ destroy() {} }),
}, { get: (t, p) => (p in t ? t[p] : p === 'then' ? undefined : generic()) });

function installDevice() {
  Object.assign(globalThis, { GPUShaderStage: { VERTEX: 1, FRAGMENT: 2, COMPUTE: 4 }, GPUBufferUsage: { MAP_READ: 1, MAP_WRITE: 2, COPY_SRC: 4, COPY_DST: 8, INDEX: 16, VERTEX: 32, UNIFORM: 64, STORAGE: 128, INDIRECT: 256, QUERY_RESOLVE: 512 }, GPUMapMode: { READ: 1, WRITE: 2 } });
  Object.defineProperty(globalThis, 'navigator', { value: { gpu: { requestAdapter: async () => ({ limits, features: new Set(), requestDevice: async () => device }) } }, configurable: true });
}

/*
 * The SHA-256 of every ocean shader module's source (the modules named oX…)
 * that createGpuModel compiles at N = 8 with `ocean` as the ocean's options
 * on top of `base`, with the engine kept in the checkout at `root`. Nothing
 * runs: a stand-in device records the sources.
 */
export async function oceanShaderDigest(ocean, { base = {}, root = HERE } = {}) {
  installDevice();
  const { createGpuModel } = await import(`${root}js/gpu/model.gpu.js`);
  const { Grid } = await import(`${root}js/grid.module.js`);
  const { sigmaInterfaces } = await import(`${root}js/dynamics/sigmaCore.module.js`);
  const { syntheticTopography } = await import(`${root}js/geography.module.js`);
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  modules = [];
  const log = console.log;
  console.log = () => {};
  try { await createGpuModel(new Grid(8), { topography, levels: sigmaInterfaces('bl36'), ocean: { everySteps: 1, ...base, ...ocean } }); } finally { console.log = log; }
  const hash = createHash('sha256');
  let count = 0;
  for (const { label, code } of modules) if (/^o[A-Z]/.test(label)) { hash.update(`${label}\n${code}\n`); count++; }
  return `${count}:${hash.digest('hex').slice(0, 16)}`;
}

const SEA_OPTIONS = [
  {},
  { minimumThickness: 10, shallowestMixedDepth: 10, refillFloor: 5, refillRate: 3 / 86400, shearEntrainment: 0.65, retreatDaily: true, retreatLight: 1, interiorShearMixing: false },
  { shearMixing: true },
];

/*
 * The SHA-256 of h, u, Q, W and η of the CPU ocean after `steps` steps (48 by default) at
 * N = 8 on a sea with a continent, under ocean option set SEA_OPTIONS[which]
 * with `ocean` on top, wind and heat flux fixed.
 */
export async function stateDigest(which, ocean = {}, { steps = 48, root = HERE } = {}) {
  const { createOcean } = await import(`${root}js/ocean/layered.module.js`);
  const { mesh, C, E, DEG, zonalWindOnEdges } = await import(`${root}test/helpers/layered.mjs`);
  const land = new Uint8Array(C), edgeOcean = new Uint8Array(E);
  for (let i = 0; i < C; i++) land[i] = Math.cos(mesh.lonCell[i]) > 0.4 && Math.abs(mesh.latCell[i]) < 55 * DEG ? 1 : 0;
  for (let e = 0; e < E; e++) edgeOcean[e] = land[mesh.cellsOnEdge[2 * e]] || land[mesh.cellsOnEdge[2 * e + 1]] ? 0 : 1;
  const geography = { land, elevation: Float64Array.from(land, (l) => (l ? 300 : -4000)), edgeOcean };
  const model = createOcean(mesh, { everySteps: 1, geography, bathymetry: new Float64Array(C).fill(4000), ...SEA_OPTIONS[which], ...ocean });
  const surfaceT = Float64Array.from(mesh.latCell, (l) => 300 - 25 * Math.sin(l) ** 2), ice = new Float64Array(C), flux = Float64Array.from(mesh.latCell, (l) => 80 * Math.cos(2 * l));
  model.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, (lat) => -0.12 * Math.cos(3 * lat));
  for (let n = 0; n < steps; n++) model.advance(Float64Array.from(surfaceT), ice, flux, wind, 2700, undefined, (n + 1) * 2700);
  const hash = createHash('sha256');
  for (const a of [model.h, model.u, model.Q, model.W, model.eta]) hash.update(Buffer.from(a.buffer, a.byteOffset, a.byteLength));
  return hash.digest('hex').slice(0, 16);
}

export const SEA_CONFIGS = SEA_OPTIONS.length;
