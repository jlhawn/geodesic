// The parts of the verda CLI's JSON (verda --agent ... -o json) that
// scripts/verdaRelaunch.sh acts on. As a command it reads that JSON on
// stdin and prints one line, or nothing when there is nothing to report:
//   instance <hostname>  from `vm list`: the id, status and OS volume
//                        ('-' when the list leaves it out) of the instance
//                        of that hostname to act on: a running one before
//                        one starting or stopping, before an offline one,
//                        before any other
//   osvolume             from `vm describe`: the instance's OS volume
//   volume <id>          from `volume list`: that volume's status
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const RANK = { running: 0, new: 1, ordered: 1, provisioning: 1, validating: 1, deleting: 1, offline: 2 };

export function listed(json) {
  if (Array.isArray(json)) return json;
  for (const key of ['instances', 'volumes', 'data', 'items']) if (Array.isArray(json?.[key])) return json[key];
  return Object.values(json ?? {}).find(Array.isArray) ?? [];
}

export function pickInstance(json, hostname) {
  const mine = listed(json).filter((instance) => instance && instance.hostname === hostname);
  return mine.sort((a, b) => (RANK[a.status] ?? 3) - (RANK[b.status] ?? 3))[0] ?? null;
}

export function osVolumeOf(instance) {
  if (!instance) return null;
  if (instance.instance && !instance.id) return osVolumeOf(instance.instance);
  if (instance.os_volume_id) return instance.os_volume_id;
  const volume = (Array.isArray(instance.volumes) ? instance.volumes : []).find((v) => v && v.is_os_volume);
  return volume ? volume.id : null;
}

export function volumeStatus(json, id) {
  const volume = listed(json).find((v) => v && v.id === id);
  return volume ? volume.status : null;
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const [mode, argument] = process.argv.slice(2);
  const text = readFileSync(0, 'utf8').trim();
  const json = text ? JSON.parse(text) : [];
  if (mode === 'instance') {
    const instance = pickInstance(json, argument);
    if (instance) console.log(`${instance.id} ${instance.status} ${osVolumeOf(instance) ?? '-'}`);
  } else if (mode === 'osvolume') {
    const volume = osVolumeOf(json);
    if (volume) console.log(volume);
  } else if (mode === 'volume') {
    const status = volumeStatus(json, argument);
    if (status) console.log(status);
  } else throw new Error(`unknown mode ${mode}: instance <hostname>, osvolume or volume <id>`);
}
