/*
 * The layered soil's WGSL (land.soil 'layers', js/physics/land.module.js):
 * the constants of both classes of column, soil and firn, and the
 * implicit conduction step that replaces the slab's `T += dt * net / LANDC`
 * in the surface kernel. The layers' changes are solved from the old
 * temperatures' heat fluxes, so rounding in float32 stays at the size of
 * the changes. The step is straight-line code on scalars, one triple per
 * layer, with no arrays, loops or function calls.
 */
export function soilConstants(column) {
  return `${[['SOIL', column.soil], ['FIRN', column.firn]].map(([tag, c]) => [
    ...c.capacity.map((value, k) => `const ${tag}_C${k + 1}: f32 = ${value};`),
    ...c.conductance.map((value, k) => `const ${tag}_G${k + 1}: f32 = ${value};`),
  ].join(' ')).join('\n')}\n`;
}

export const SOIL_SKIN = 'colC1';

export function soilStep(column) {
  const n = column.depths.length, pick = (tag, k) => `select(SOIL_${tag}${k}, FIRN_${tag}${k}, onIceSheet)`;
  const lines = [`let ${SOIL_SKIN} = ${pick('C', 1)};`];
  if (n === 1) {
    lines.push(`T += net / (${SOIL_SKIN} / dt);`);
    return lines.join('\n    ');
  }
  const k1 = Array.from({ length: n - 1 }, (_, k) => k + 1);
  lines.push('{');
  for (const k of k1) lines.push(`  let colG${k} = ${pick('G', k)};`);
  for (let k = 2; k <= n; k++) lines.push(`  let colT${k} = PH[PH_SOILT${k} + i];`);
  lines.push(`  let colF1 = colG1 * (T - colT2);`);
  for (let k = 2; k < n; k++) lines.push(`  let colF${k} = colG${k} * (colT${k} - colT${k + 1});`);
  lines.push(`  let colB1 = ${SOIL_SKIN} / dt + colG1; let colS1 = (net - colF1) / colB1; let colP1 = colG1 / colB1;`);
  for (let k = 2; k < n; k++) {
    lines.push(`  let colB${k} = ${pick('C', k)} / dt + colG${k - 1} + colG${k} - colG${k - 1} * colP${k - 1}; let colS${k} = (colF${k - 1} - colF${k} + colG${k - 1} * colS${k - 1}) / colB${k}; let colP${k} = colG${k} / colB${k};`);
  }
  lines.push(`  let colB${n} = ${pick('C', n)} / dt + colG${n - 1} - colG${n - 1} * colP${n - 1}; let colE${n} = (colF${n - 1} + colG${n - 1} * colS${n - 1}) / colB${n};`);
  for (let k = n - 1; k >= 1; k--) lines.push(`  let colE${k} = colS${k} + colP${k} * colE${k + 1};`);
  lines.push('  T += colE1;');
  for (let k = 2; k <= n; k++) lines.push(`  PH[PH_SOILT${k} + i] = colT${k} + colE${k};`);
  lines.push('}');
  return lines.join('\n    ');
}
