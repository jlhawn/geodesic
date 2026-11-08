import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, dirname, resolve, relative } from 'node:path';
import { spawnSync } from 'node:child_process';
import { pathToFileURL } from 'node:url';
import { Grid } from '../js/grid.module.js';
import { gridHeader } from '../js/stateGrid.module.js';
import { FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { encodeState } from '../js/stateFile.module.js';
import { unreadOptions, withoutOptionGroups } from '../scripts/optionEnvironment.mjs';
import { headerOptions } from '../scripts/stateOptions.mjs';

const scripts = new URL('../scripts/', import.meta.url).pathname, dir = mkdtempSync(join(tmpdir(), 'option-environment-'));
after(() => rmSync(dir, { recursive: true, force: true }));
const run = (file, args, env) => spawnSync(process.execPath, [file, ...args], { env: { ...withoutOptionGroups(), ...env }, encoding: 'utf8' });

const files = readdirSync(scripts, { recursive: true }).filter((name) => name.endsWith('.mjs')).map((name) => join(scripts, name));
const text = new Map(files.map((file) => [file, readFileSync(file, 'utf8')]));
const imports = (file) => [...text.get(file).matchAll(/\bfrom\s*['"](\.{1,2}\/[^'"]+\.mjs)['"]/g)].map((match) => resolve(dirname(file), match[1])).filter((path) => text.has(path));
const script = (file) => /process\.argv/.test(text.get(file));
const readsGroups = (file, seen = new Set()) => !seen.has(file) && (seen.add(file), /JSON\.parse\(\s*(process\.env|env\()/.test(text.get(file)) || imports(file).some((helper) => !script(helper) && readsGroups(helper, seen)));
const imported = new Set(files.flatMap(imports));
const groupsRead = (file, seen = new Set()) => {
  if (seen.has(file)) return [];
  seen.add(file);
  const own = [...text.get(file).matchAll(/JSON\.parse\(\s*process\.env\.([A-Z][A-Z0-9_]*)|\benv\(\s*'([A-Z][A-Z0-9_]*)'|\[\s*'[a-z][A-Za-z]*'\s*,\s*'([A-Z][A-Z0-9_]*)'\s*\]/g)].map((match) => match[1] ?? match[2] ?? match[3]);
  return [...new Set([...own, ...imports(file).filter((helper) => !script(helper)).flatMap((helper) => groupsRead(helper, seen))])];
};

test('every script that reads option groups from the environment, itself or through a helper it imports, refuses one under a name it does not read before it does anything else, and names only that one beside every group it reads', () => {
  const entries = files.filter((file) => (script(file) || !imported.has(file)) && readsGroups(file)).map((file) => relative(scripts, file));
  for (const name of ['spinup.mjs', 'fixedSst.mjs', 'oceanSpinup.mjs', 'upperAtmosphere.mjs', 'figures/stateMaps.mjs']) assert.ok(entries.includes(name), name);
  for (const [name, group] of [['spinup.mjs', 'DYNAMICS'], ['fixedSst.mjs', 'GIVEN'], ['fixedSst.mjs', 'MOIST'], ['figures/stateMaps.mjs', 'BOUNDARY_LAYER'], ['armMetrics.mjs', 'SURFACE']]) assert.ok(groupsRead(join(scripts, name)).includes(group), `${name} ${group}`);
  const left = entries.filter((name) => {
    const result = run(join(scripts, name), ['/nonexistent/state.bin'], { ...Object.fromEntries(groupsRead(join(scripts, name)).map((group) => [group, '{}'])), UNREAD_GROUP: '{"x":1}' });
    return result.status === 0 || !result.stderr.includes(`${name.split('/').pop()} does not read UNREAD_GROUP, given as an option group`);
  });
  assert.deepEqual(left, []);
});

test('a script reads the names its text shows and those of the helpers it imports that are not scripts of their own, and refuses an option group under any other name only when node runs it', () => {
  const check = pathToFileURL(join(scripts, 'optionEnvironment.mjs')).href;
  writeFileSync(join(dir, 'helper.mjs'), "export const helped = () => JSON.parse(process.env.HELPED ?? '{}');\n");
  writeFileSync(join(dir, 'tool.mjs'), "export const tooled = 1;\nif (process.argv[2] === 'tool') console.log(JSON.parse(process.env.TOOLED ?? '{}'));\n");
  writeFileSync(join(dir, 'main.mjs'), `import { refuseUnreadOptions } from '${check}';\nimport { helped } from './helper.mjs';\nimport { tooled } from './tool.mjs';\n\nrefuseUnreadOptions(import.meta.url);\nconst env = (name) => JSON.parse(process.env[name] ?? '{}');\nconsole.log('ran', JSON.stringify([process.env.OWN ?? null, env('QUOTED'), helped(), tooled]));\n`);
  writeFileSync(join(dir, 'outer.mjs'), "import './main.mjs';\n");
  const main = pathToFileURL(join(dir, 'main.mjs'));
  assert.deepEqual(unreadOptions(main, { OWN: '{}', QUOTED: '{"a":1}', HELPED: '{}', TOOLED: '{}', OTHER: ' {"b":2}', PLAIN: 'x', LIST: '[1]' }), ['OTHER', 'TOOLED']);
  const given = { OWN: '{}', QUOTED: '{"a":1}', HELPED: '{"c":3}' }, ran = run(main.pathname, [], given);
  assert.equal(ran.status, 0, ran.stderr);
  assert.match(ran.stdout, /^ran \["\{\}",\{"a":1\},\{"c":3\},1\]/);
  const refused = run(main.pathname, [], { ...given, TOOLED: '{}' });
  assert.notEqual(refused.status, 0);
  assert.match(refused.stderr, /main\.mjs does not read TOOLED, given as an option group: unset or rename it/);
  assert.doesNotMatch(refused.stdout, /ran/);
  const outer = run(join(dir, 'outer.mjs'), [], { TOOLED: '{}' });
  assert.equal(outer.status, 0, outer.stderr);
});

test('fixedSst.mjs refuses a state whose header records no given options, and an ocean group, which a run over a fixed sea does not read, before it builds anything, and takes GIVEN in place of the record with each group the environment gives over it, but not beside a record it differs from', () => {
  const mesh = new Grid(16), file = join(dir, 'headerOnly_day0000.bin');
  writeFileSync(file, encodeState({ N: 16, ...gridHeader(mesh), ...FRESH_GEOGRAPHY, terrain: true, K: 48, day: 0, time: 0, options: { radiation: {}, moist: {}, land: {} } }));
  const fixed = (env, state = file) => {
    const environment = { ...withoutOptionGroups(), ...env };
    delete environment.GPU_QUEUE_HELD;
    return spawnSync(process.execPath, [join(scripts, 'fixedSst.mjs'), state, dir, 'arm'], { env: environment, encoding: 'utf8' });
  };
  const plain = fixed({});
  assert.notEqual(plain.status, 0);
  assert.match(plain.stderr, /records no given options .* GIVEN \(JSON, \{\} for today's defaults\) stands in for them/);
  assert.doesNotMatch(plain.stdout, /options:/);
  const ocean = fixed({ GIVEN: '{}', OCEAN: '{"mixing":"kpp"}' });
  assert.notEqual(ocean.status, 0);
  assert.match(ocean.stderr, /fixedSst\.mjs does not read OCEAN, given as an option group/);
  const standIn = fixed({ GIVEN: JSON.stringify({ moist: { pcapeMinimum: 0.05 }, orography: false, ocean: { mixing: 'kpp' } }), MOIST: '{"bergeronAutoconversion":false}', OROGRAPHY: '{"constants":"cy41r2"}' });
  const line = standIn.stdout.split('\n').find((text) => text.startsWith('options: '));
  assert.ok(line, standIn.stdout + standIn.stderr);
  assert.match(line, /^options: GIVEN, the header recording none, MOIST, OROGRAPHY over it: /);
  assert.deepEqual(JSON.parse(line.slice(line.indexOf('{'))), { moist: { pcapeMinimum: 0.05, bergeronAutoconversion: false }, orography: { constants: 'cy41r2' } });
  const recorded = join(dir, 'recorded_day0000.bin');
  writeFileSync(recorded, encodeState({ N: 16, ...gridHeader(mesh), ...FRESH_GEOGRAPHY, terrain: true, K: 48, day: 0, time: 0, ...headerOptions({ moist: { pcapeMinimum: 0 } }) }));
  const differing = fixed({ GIVEN: '{"moist":{"pcapeMinimum":0.05}}' }, recorded);
  assert.notEqual(differing.status, 0);
  assert.match(differing.stderr, /this header records them, and GIVEN differs from them in radiation, moist, land: unset GIVEN/);
  assert.doesNotMatch(differing.stdout, /options:/);
});
