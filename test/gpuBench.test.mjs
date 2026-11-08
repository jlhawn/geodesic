import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, readFileSync, existsSync, mkdtempSync, mkdirSync, rmSync, chmodSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { withoutOptionGroups } from '../scripts/optionEnvironment.mjs';

const bench = new URL('../scripts/gpuBench.sh', import.meta.url).pathname, dir = mkdtempSync(join(tmpdir(), 'gpu-bench-'));
after(() => rmSync(dir, { recursive: true, force: true }));
const bin = join(dir, 'bin'), trace = join(dir, 'trace.txt');
mkdirSync(bin);
writeFileSync(join(bin, 'node'), `#!/bin/bash
case "$1" in
  scripts/spinup.mjs|scripts/profileGpu.mjs|scripts/paceGpu.mjs) ;;
  *) exec '${process.execPath}' "$@";;
esac
echo "$1 $(env | sed -n 's/^\\([A-Za-z_][A-Za-z0-9_]*\\)=\\([[:space:]]*{.*\\)/\\1=\\2/p' | sort | tr '\\n' '\\t')" >> '${trace}'
[ "$1" = scripts/spinup.mjs ] && { : > "$OUT/\${TAG}_day0000_step\${STOP_AFTER_STEPS}.bin"; echo 'after 1 s of setup' > "$OUT/$TAG.log"; }
echo "$1 ran"
`);
chmodSync(join(bin, 'node'), 0o755);
const run = (env) => {
  writeFileSync(trace, '');
  return spawnSync('bash', [bench], { env: { ...withoutOptionGroups(), PATH: `${bin}:${process.env.PATH}`, NS: '', FRESH: '16', FRESH_STEPS: '2', ...env }, encoding: 'utf8' });
};
const stages = () => readFileSync(trace, 'utf8').split('\n').filter(Boolean).map((line) => {
  const at = line.indexOf(' ');
  return [line.slice(0, at), Object.fromEntries(line.slice(at + 1).split('\t').filter(Boolean).map((entry) => [entry.slice(0, entry.indexOf('=')), entry.slice(entry.indexOf('=') + 1)]))];
});

test("gpuBench.sh gives each stage only the option groups its script reads, a value's leading blank kept, takes FRESH_OCEAN as the fresh start's OCEAN, and names in its report the groups a script runs without", () => {
  const out = join(dir, 'plumbed'), ocean = '{"everySteps":2}', moist = ' {"pcapeMinimum":0}';
  const result = run({ OUT: out, FRESH_OCEAN: ocean, DYNAMICS: '{"thetaFace":"arakawaSuarez"}', MOIST: moist });
  assert.equal(result.status, 0, result.stderr);
  const ran = stages();
  assert.deepEqual(ran.map(([script]) => script), ['scripts/spinup.mjs', 'scripts/profileGpu.mjs', 'scripts/profileGpu.mjs', 'scripts/paceGpu.mjs', 'scripts/paceGpu.mjs']);
  assert.deepEqual(ran[0][1], { DYNAMICS: '{"thetaFace":"arakawaSuarez"}', MOIST: moist, OCEAN: ocean });
  for (const [script, groups] of ran.slice(1)) assert.deepEqual(groups, { MOIST: moist, OCEAN: ocean }, script);
  const report = readFileSync(join(out, 'report.txt'), 'utf8');
  for (const script of ['profileGpu', 'paceGpu']) assert.match(report, new RegExp(`^${script}\\.mjs runs without the option groups it does not read: DYNAMICS$`, 'm'));
  assert.doesNotMatch(report, /spinup\.mjs runs without|FAILED/);
});

test('gpuBench.sh refuses an option group no stage reads before it runs or writes anything', () => {
  const out = join(dir, 'refused'), result = run({ OUT: out, DYNAMIC: '{"thetaFace":"arakawaSuarez"}', MOIST: '{}' });
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /no stage reads these option groups, so the bench would run without them: DYNAMIC; unset or rename them/);
  assert.equal(existsSync(out), false);
  assert.deepEqual(stages(), []);
});
