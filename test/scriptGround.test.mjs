import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';

const scripts = new URL('../scripts/', import.meta.url);
const BUILDERS = ['createModel', 'createGpuModel', 'createParallelModel', 'createGeography'];
const GROUND = ['stateGeography', 'savedGeography'];
const NAME = /(?<![\w$])(?:(?<=\.\.\.)|(?<!\.))[A-Za-z_$][\w$]*/g;

function closing(text, i) {
  for (let depth = 1; depth && i < text.length; i++) {
    if ('([{'.includes(text[i])) depth++;
    else if (')]}'.includes(text[i])) depth--;
  }
  return i;
}

function parts(text) {
  const out = [];
  let depth = 0, from = 0;
  for (let i = 0; i <= text.length; i++) {
    const c = text[i];
    if ('([{'.includes(c)) depth++;
    else if (')]}'.includes(c)) depth--;
    else if (i === text.length || (c === ',' && depth === 0)) {
      if (text.slice(from, i).trim()) out.push(text.slice(from, i).trim());
      from = i + 1;
    }
  }
  return out;
}

function expression(text, i) {
  let depth = 0, j = i;
  for (; j < text.length; j++) {
    const c = text[j];
    if ('([{'.includes(c)) depth++;
    else if (')]}'.includes(c)) { if (!depth) break; depth--; }
    else if ((c === ',' || c === ';') && !depth) break;
  }
  return text.slice(i, j).trim();
}

function definitions(text) {
  const defs = new Map(), add = (name, value) => defs.set(name, [...(defs.get(name) ?? []), value]);
  for (const match of text.matchAll(/(?<![\w$.])([A-Za-z_$][\w$]*)(?:\.[\w$]+)*\s*=(?![=>])/g)) add(match[1], expression(text, match.index + match[0].length));
  for (const match of text.matchAll(/[{[]([^{}[\]]*)[}\]]\s*=(?![=>])/g)) {
    const value = expression(text, match.index + match[0].length);
    for (const part of parts(match[1])) {
      const name = /([A-Za-z_$][\w$]*)\s*(?:=[\s\S]*)?$/.exec(part);
      if (name) add(name[1], value);
    }
  }
  for (const match of text.matchAll(/\bfunction\s*\*?\s*([A-Za-z_$][\w$]*)\s*\(/g)) {
    const body = text.indexOf('{', closing(text, match.index + match[0].length));
    add(match[1], text.slice(match.index, closing(text, body + 1)));
  }
  return defs;
}

function grounded(value, defs, ground, seen = new Set()) {
  value = value.trim();
  if (value.startsWith('{') && closing(value, 1) === value.length) {
    const members = parts(value.slice(1, -1));
    return members.length > 0 && members.every((part) => grounded(part.startsWith('...') ? part.slice(3) : part.replace(/^[A-Za-z_$][\w$]*\s*:/, ''), defs, ground, new Set(seen)));
  }
  if (new RegExp(`\\b(${ground.join('|')})\\(`).test(value)) return true;
  for (const { 0: name, index } of value.matchAll(NAME)) {
    if (seen.has(name) || (/[{,]\s*$/.test(value.slice(0, index)) && /^\s*:/.test(value.slice(index + name.length)))) continue;
    seen.add(name);
    if ((defs.get(name) ?? []).some((def) => grounded(def, defs, ground, seen))) return true;
  }
  return false;
}

function entries(value, defs, seen = new Set()) {
  value = value.trim();
  if (value.startsWith('{') && closing(value, 1) === value.length) {
    return parts(value.slice(1, -1)).flatMap((part) => {
      if (part.startsWith('...')) return entries(part.slice(3), defs, seen);
      const entry = /^([A-Za-z_$][\w$]*)\s*(?::([\s\S]*))?$/.exec(part);
      return entry ? [{ key: entry[1], value: (entry[2] ?? entry[1]).trim() }] : [];
    });
  }
  if (!/^[A-Za-z_$][\w$]*$/.test(value) || seen.has(value)) return [];
  seen.add(value);
  return (defs.get(value) ?? []).flatMap((def) => entries(def, defs, seen));
}

function calls(text, names) {
  const found = [];
  for (const match of text.matchAll(new RegExp(`\\b(${names.join('|')})\\(`, 'g'))) {
    if (/^\s*(\/\/|\/?\*)/.test(text.slice(text.lastIndexOf('\n', match.index) + 1, match.index))) continue;
    const start = match.index + match[0].length, end = closing(text, start);
    found.push({ name: match[1], line: text.slice(0, match.index).split('\n').length, text: text.slice(match.index, end), args: parts(text.slice(start, end - 1)) });
  }
  return found;
}

function misplaced(file, text) {
  const left = [], defs = definitions(text), ground = [...GROUND], builders = new Map(BUILDERS.map((name) => [name, name]));
  for (const match of text.matchAll(new RegExp(`\\b(${BUILDERS.join('|')})\\s*(?:\\bas\\b|:)\\s*([A-Za-z_$][\\w$]*)`, 'g'))) if (!BUILDERS.includes(match[2])) builders.set(match[2], match[1]);
  for (const match of text.matchAll(new RegExp(`([A-Za-z_$][\\w$]*)\\s*:\\s*(${BUILDERS.join('|')})\\b(?!\\s*\\()`, 'g'))) if (!BUILDERS.includes(match[1])) builders.set(match[1], match[2]);
  for (const [alias, values] of defs) for (const value of values) if (BUILDERS.includes(value)) builders.set(alias, value);
  for (const match of text.matchAll(/\bimport\s*\{([^}]*)\}\s*from\s*['"](\.{1,2}\/[^'"]+\.mjs)['"]/g)) {
    const module = new URL(match[2], new URL(file, scripts)), theirs = definitions(readFileSync(module, 'utf8'));
    for (const [, name, local] of parts(match[1]).map((part) => /^([\w$]+)(?:\s+as\s+([\w$]+))?$/.exec(part)).filter(Boolean)) {
      if ((theirs.get(name) ?? []).some((def) => grounded(def, theirs, GROUND))) ground.push(local ?? name);
    }
  }
  if (/topography_0p25\.bin|topographyFromInt16\(/.test(text)) left.push(`${file} reads the raster's topography itself`);
  for (const call of calls(text, [...builders.keys()])) {
    const where = `${file}:${call.line} ${call.name}`;
    if (builders.get(call.name) === 'createGeography') {
      if (call.args.length < 3) left.push(`${where} without the state's geography options`);
      for (const [i, role] of [[1, 'topography'], [2, 'geography']]) if (call.args[i] !== undefined && !grounded(call.args[i], defs, ground)) left.push(`${where}'s ${role} ${call.args[i]} is not the state's own`);
      continue;
    }
    const named = entries(call.args[1] ?? '', defs).filter((entry) => entry.key === 'topography' || entry.key === 'geography');
    if ((named.length || /\btopography\b/.test(call.text)) && !named.some((entry) => entry.key === 'geography')) left.push(`${where} without the state's geography`);
    for (const entry of named) if (!grounded(entry.value, defs, ground)) left.push(`${where}'s ${entry.key} ${entry.value} is not the state's own`);
  }
  return left;
}

test('every script that reads a saved state builds its models and geographies on the topography and geography that stateGeography or savedGeography gives for that state, never on the raster of data/topography_0p25.bin, as far as the scripts\' text shows: options that resolve to no object literal in the file, such as a call\'s result, are taken on trust, a model given no topography is an aquaplanet, and a value is the state\'s own when any name in it leads there through any assignment in its file', () => {
  const left = [];
  for (const file of readdirSync(scripts, { recursive: true }).filter((name) => name.endsWith('.mjs')).sort()) {
    const text = readFileSync(new URL(file, scripts), 'utf8');
    if (/\b(decodeState|readState|stateGrid)\(/.test(text)) left.push(...misplaced(file, text));
  }
  assert.deepEqual(left, []);
});
