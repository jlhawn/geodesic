import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, readFileSync, mkdtempSync, mkdirSync, rmSync, chmodSync, existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { listed, pickInstance, osVolumeOf, volumeStatus } from '../scripts/verdaInstances.mjs';

const root = new URL('..', import.meta.url).pathname;
let dir;
before(() => { dir = mkdtempSync(join(tmpdir(), 'verdaRelaunch-')); });
after(() => rmSync(dir, { recursive: true, force: true }));

const FAKE = `#!/bin/bash
echo "$*" >> "$SCENARIO/calls"
n=$(cat "$SCENARIO/round" 2>/dev/null || echo 0)
case "$*" in *"vm list"*) n=$((n + 1)); echo $n > "$SCENARIO/round" ;; esac
answer() { if [ -f "$SCENARIO/$1.$n.json" ]; then cat "$SCENARIO/$1.$n.json"; elif [ -f "$SCENARIO/$1.json" ]; then cat "$SCENARIO/$1.json"; else echo '[]'; fi; }
case "$*" in
  *"vm list"*) answer list ;;
  *"vm describe"*) answer describe ;;
  *"volume list"*) answer volumes ;;
  *"vm create"*) echo '{"id":"i-new","status":"provisioning"}' ;;
  *"vm start"*) echo '{"status":"ok"}' ;;
esac
`;

function scenario(name, files) {
  const at = join(dir, name);
  mkdirSync(at);
  for (const [file, json] of Object.entries(files)) writeFileSync(join(at, file), JSON.stringify(json));
  writeFileSync(join(at, 'verda'), FAKE);
  chmodSync(join(at, 'verda'), 0o755);
  return at;
}
function relaunch(at, rounds, env = {}) {
  const result = spawnSync('bash', [join(root, 'scripts/verdaRelaunch.sh')], {
    cwd: root, encoding: 'utf8',
    env: { ...process.env, SCENARIO: at, VERDA: join(at, 'verda'), NAME: 'gcm64', INSTANCE_TYPE: '1A100.22V', OS: 'ubuntu-24.04-cuda-12.8-open-docker', SSH_KEY: 'key-1', STARTUP_SCRIPT: 'script-1', POLL: '0', ROUNDS: String(rounds), STATE_FILE: join(at, 'state'), LOG: join(at, 'log'), STOP_FILE: join(at, 'STOP'), ...env },
  });
  assert.equal(result.status, 0, result.stderr);
  return { calls: readFileSync(join(at, 'calls'), 'utf8').trim().split('\n'), log: readFileSync(join(at, 'log'), 'utf8') };
}
const creates = (calls) => calls.filter((call) => call.includes('vm create'));

test('the JSON helper finds the instance to act on and its OS volume in the shapes the CLI may print', () => {
  const instances = [
    { id: 'i-old', hostname: 'gcm64', status: 'discontinued' },
    { id: 'i-other', hostname: 'other', status: 'running', os_volume_id: 'vol-x' },
    { id: 'i-now', hostname: 'gcm64', status: 'running', os_volume_id: 'vol-1' },
  ];
  assert.deepEqual(listed({ instances }), instances);
  assert.deepEqual(listed({ data: instances }), instances);
  assert.deepEqual(listed(null), []);
  assert.equal(pickInstance(instances, 'gcm64').id, 'i-now');
  assert.equal(pickInstance([{ id: 'a', hostname: 'gcm64', status: 'offline' }, { id: 'b', hostname: 'gcm64', status: 'provisioning' }], 'gcm64').id, 'b');
  assert.equal(pickInstance(instances, 'absent'), null);
  assert.equal(osVolumeOf({ instance: { id: 'i', volumes: [{ id: 'vol-data', is_os_volume: false }, { id: 'vol-os', is_os_volume: true }] } }), 'vol-os');
  assert.equal(osVolumeOf({ id: 'i' }), null);
  assert.equal(volumeStatus({ volumes: [{ id: 'vol-1', status: 'detached' }] }, 'vol-1'), 'detached');
  assert.equal(volumeStatus([], 'vol-1'), null);
});

test('the relauncher remembers the OS volume, waits for it to detach after an eviction, recreates the spot instance on it and starts an offline one', () => {
  const at = scenario('evicted', {
    'list.1.json': [{ id: 'i-1', hostname: 'gcm64', status: 'running' }],
    'describe.json': { id: 'i-1', hostname: 'gcm64', status: 'running', volumes: [{ id: 'vol-1', is_os_volume: true }] },
    'list.2.json': [],
    'volumes.2.json': [{ id: 'vol-1', status: 'attached' }],
    'list.3.json': { instances: [{ id: 'i-1', hostname: 'gcm64', status: 'discontinued' }] },
    'volumes.3.json': [{ id: 'vol-1', status: 'detached' }],
    'list.4.json': [{ id: 'i-2', hostname: 'gcm64', status: 'provisioning' }],
    'list.5.json': [{ id: 'i-2', hostname: 'gcm64', status: 'offline' }],
    'list.6.json': [{ id: 'i-2', hostname: 'gcm64', status: 'running', os_volume_id: 'vol-1' }],
  });
  const { calls, log } = relaunch(at, 6);
  assert.equal(readFileSync(join(at, 'state'), 'utf8').trim(), 'vol-1');
  assert.equal(creates(calls).length, 1);
  const create = creates(calls)[0];
  for (const flag of ['--agent', 'vm create', '--kind gpu', '--instance-type 1A100.22V', '--location FIN-01', '--is-spot', '--os vol-1', '--os-volume-size 100', '--os-volume-on-spot-discontinue keep_detached', '--ssh-key key-1', '--hostname gcm64', '--startup-script script-1', '--wait', '-o json']) assert.ok(create.includes(flag), `${flag} in: ${create}`);
  assert.ok(calls.includes('--agent vm start i-2 -o json'), calls.join('\n'));
  assert.match(log, /\d{4}-\d\d-\d\d \d\d:\d\d:\d\d gcm64 is running \(i-1\)/);
  assert.match(log, /remembering OS volume vol-1 of gcm64/);
  assert.match(log, /waiting for OS volume vol-1 to detach \(it is attached\)/);
  assert.match(log, /creating gcm64: 1A100\.22V spot in FIN-01 on vol-1/);
  assert.match(log, /starting gcm64 \(i-2\)/);
  assert.equal(log.match(/remembering/g).length, 1, 'the same volume is remembered once');
});

test('with no OS volume known the first instance comes from the image, and a volume that vanished stops the relauncher from creating on another', () => {
  const first = relaunch(scenario('fresh', { 'list.json': [] }), 1);
  assert.equal(creates(first.calls).length, 1);
  assert.ok(creates(first.calls)[0].includes('--os ubuntu-24.04-cuda-12.8-open-docker'));

  const lost = scenario('lost', { 'list.json': [], 'volumes.json': [{ id: 'vol-other', status: 'detached' }] });
  writeFileSync(join(lost, 'state'), 'vol-1\n');
  const { calls, log } = relaunch(lost, 2);
  assert.equal(creates(calls).length, 0);
  assert.match(log, /OS volume vol-1 is not listed; not creating gcm64 on another/);

  const given = scenario('given', { 'list.json': [], 'volumes.json': [{ id: 'vol-9', status: 'detached' }] });
  assert.ok(creates(relaunch(given, 1, { OS_VOLUME: 'vol-9' }).calls)[0].includes('--os vol-9'));
});

test('the relauncher stops at its stop file', () => {
  const at = scenario('stopped', { 'list.json': [{ id: 'i-1', hostname: 'gcm64', status: 'running', os_volume_id: 'vol-1' }] });
  writeFileSync(join(at, 'STOP'), '');
  writeFileSync(join(at, 'calls'), '');
  const result = spawnSync('bash', [join(root, 'scripts/verdaRelaunch.sh')], { cwd: root, encoding: 'utf8', env: { ...process.env, SCENARIO: at, VERDA: join(at, 'verda'), NAME: 'gcm64', INSTANCE_TYPE: 't', SSH_KEY: 'k', STARTUP_SCRIPT: 's', POLL: '0', STATE_FILE: join(at, 'state'), LOG: join(at, 'log'), STOP_FILE: join(at, 'STOP') } });
  assert.equal(result.status, 0);
  assert.equal(readFileSync(join(at, 'calls'), 'utf8'), '');
  assert.match(readFileSync(join(at, 'log'), 'utf8'), /stopped at .*STOP/);
  assert.ok(!existsSync(join(at, 'state')));
});
