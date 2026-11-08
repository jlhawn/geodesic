import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';

const scripts = new URL('../scripts/', import.meta.url);

function calls(text, names) {
  const found = [];
  for (const match of text.matchAll(new RegExp(`\\b(${names.join('|')})\\(`, 'g'))) {
    let depth = 1, i = match.index + match[0].length, args = 1;
    while (depth && i < text.length) {
      const c = text[i++];
      if ('([{'.includes(c)) depth++;
      else if (')]}'.includes(c)) depth--;
      else if (c === ',' && depth === 1) args++;
    }
    found.push({ name: match[1], line: text.slice(0, match.index).split('\n').length, text: text.slice(match.index, i), args });
  }
  return found;
}

test('every script that finds a saved state\'s grid builds its models and geography on that state\'s own ground, never on the raster of data/topography_0p25.bin', () => {
  const left = [];
  for (const name of readdirSync(scripts, { recursive: true }).filter((file) => file.endsWith('.mjs')).sort()) {
    const text = readFileSync(new URL(name, scripts), 'utf8');
    if (!/\bstateGrid\(/.test(text)) continue;
    if (/topography_0p25\.bin|topographyFromInt16\(/.test(text)) left.push(`${name} reads the raster's topography itself`);
    for (const call of calls(text, ['createModel', 'createGpuModel', 'createParallelModel'])) {
      if (/\btopography\b/.test(call.text) && !/\bgeography\b/.test(call.text)) left.push(`${name}:${call.line} ${call.name} without the state's geography`);
    }
    for (const call of calls(text, ['createGeography'])) if (call.args < 3) left.push(`${name}:${call.line} createGeography without the state's geography options`);
  }
  assert.deepEqual(left, []);
});
