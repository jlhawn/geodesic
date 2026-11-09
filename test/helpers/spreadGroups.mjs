/*
 * S9's groups for the distribution checks: every cell, or edge, in exactly
 * one surface (sea, land, or for an edge between the two, coast), and a 3-D
 * field split by its layers' reference pressure, σ·1000 hPa, at 700 and
 * 100 hPa. A group's gap at a cell is the signed gap of largest magnitude
 * over the group's layers, and its reference is the CPU field over the same
 * cells and layers. bandLayers puts every layer in one band and throws on a
 * sigma that is not finite; spreadGroups throws on bands that miss or repeat
 * a layer and on a field that is not their layers times the cells or edges.
 * groupChecks requires every tier of S9 in each group's bounds. The tests'
 * bounds were built from CPU twins grouped by these same functions, so a
 * change here needs the bounds rebuilt.
 */
import { hemispheres, columnLargest, magnitude, distributionCheck, S9_TIERS } from './distribution.mjs';

const SURFACES = ['sea', 'land', 'coast'];
export const BANDS = [['low', (sigma) => sigma >= 0.7], ['mid', (sigma) => sigma >= 0.1 && sigma < 0.7], ['high', (sigma) => sigma < 0.1]];

export function groupLayout(mesh, { land = null, edges = false } = {}) {
  const n = edges ? mesh.nEdges : mesh.nCells, wet = (i) => (land?.[i] ? 1 : 0);
  const kinds = Uint8Array.from({ length: n }, (_, x) => { if (!edges) return wet(x); const a = wet(mesh.cellsOnEdge[2 * x]), b = wet(mesh.cellsOnEdge[2 * x + 1]); return a === b ? a : 2; });
  const sides = hemispheres(edges ? { xCell: mesh.xEdge, nCells: n } : mesh), members = {};
  kinds.forEach((kind, x) => (members[SURFACES[kind]] ??= []).push(x));
  return { n, sides, members };
}

export function bandLayers(sigmaMid, layers) {
  if (layers === 1) return { all: [0] };
  if (!Number.isInteger(layers) || layers < 1) throw new Error(`${layers} layers`);
  const sigma = Array.from({ length: layers }, (_, k) => sigmaMid?.[k]), bad = sigma.findIndex((s) => !Number.isFinite(s));
  if (bad >= 0) throw new Error(`layer ${bad} of ${layers} has a sigma of ${sigma[bad]}`);
  return Object.fromEntries(BANDS.map(([name, inBand]) => [name, sigma.flatMap((s, k) => (inBand(s) ? [k] : []))]).filter(([, ks]) => ks.length));
}

export function spreadGroups(cpu, gpu, layout, bands) {
  const { n, sides, members } = layout, held = Object.values(bands).flat().sort((a, b) => a - b), layers = held.length;
  if (!layers || held.some((k, j) => k !== j)) throw new Error(`the bands hold layers ${held.join(', ')}, not each of 0 to ${layers - 1} once`);
  if (cpu.length !== layers * n || gpu.length !== layers * n) throw new Error(`fields of ${cpu.length} and ${gpu.length} values for ${layers} layers of ${n}`);
  const gaps = Float64Array.from({ length: layers * n }, (_, x) => gpu[x] - cpu[x]), groups = {};
  for (const [band, ks] of Object.entries(bands)) {
    const largest = columnLargest(gaps, n, ks);
    for (const [surface, cells] of Object.entries(members)) {
      groups[band === 'all' ? surface : `${surface} ${band}`] = { cells, gaps: Float64Array.from(cells, (x) => largest[x]), where: Int8Array.from(cells, (x) => sides[x]), reference: magnitude(cells.flatMap((x) => ks.map((k) => cpu[k * n + x]))) };
    }
  }
  return groups;
}

export function groupChecks(groups, spreads, { require = S9_TIERS } = {}) {
  const failed = [], rows = [];
  for (const [name, group] of Object.entries(groups)) {
    if (!spreads[name]) { failed.push(`${name}: no bounds for this group`); continue; }
    const check = distributionCheck(group.gaps, group.where, spreads[name], group.reference, { cells: group.cells, require });
    failed.push(...check.failed.map((text) => `${name}: ${text}`));
    rows.push(`${name}: ${check.text}`);
  }
  for (const name of Object.keys(spreads)) if (!groups[name]) failed.push(`${name}: bounds for a group the run does not have`);
  return { failed, rows };
}
