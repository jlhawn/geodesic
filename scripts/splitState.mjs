// Usage: node scripts/splitState.mjs runs/<name>.json|.bin[.gz] [chunkMiB]
//        node scripts/splitState.mjs runs/<name>.parts.json
// Writes <name>.json.gz.000 (or .bin.gz.000), .001, … of at most chunkMiB
// each and the <name>.parts.json manifest the page loads them from, each
// part listed with its SHA-256 and the run with its id, the SHA-256 of
// those hashes one per line. A binary state goes in without its SST ring
// (js/sstMeans.module.js), which the page keeps no recorder for, and with
// RADIATION and MOIST (JSON, as scripts/spinup.mjs takes them) recorded in
// its header as the options the page's model starts it with. Given a
// manifest, rewrites it with the hashes of the parts it lists.
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';
import { decodeState, encodeState, f64Names } from '../js/stateFile.module.js';
import { SST_FIELDS, withoutSstRing } from '../js/sstMeans.module.js';

const [file, mib = '10'] = process.argv.slice(2);
if (!file) { console.error('usage: node scripts/splitState.mjs <state.json|state.bin[.gz]|name.parts.json> [chunkMiB]'); process.exit(1); }
const sha256 = (bytes) => createHash('sha256').update(bytes).digest('hex');
const manifest = (name, description, parts) => JSON.stringify({ name, ...(description ? { description } : {}), id: sha256(parts.map((part) => part.sha256).join('\n')), parts }, null, 1) + '\n';

if (file.endsWith('.parts.json')) {
  const listed = JSON.parse(readFileSync(file, 'utf8')), dir = file.replace(/[^/]*$/, '');
  const parts = listed.parts.map((part) => { const bytes = readFileSync(dir + part.file); return { file: part.file, bytes: bytes.length, sha256: sha256(bytes) }; });
  writeFileSync(file, manifest(listed.name ?? file.replace(/.*\//, '').replace(/\.parts\.json$/, ''), listed.description, parts));
  console.log(`${file}: ${parts.length} parts hashed`);
} else {
  const chunk = Number(mib) * 1048576;
  const kind = /\.bin(?:\.gz)?$/.test(file) ? 'bin' : 'json';
  const base = file.replace(/\.(?:json|bin)(?:\.gz)?$/, '');
  const options = Object.fromEntries([['radiation', process.env.RADIATION], ['moist', process.env.MOIST]].filter(([, json]) => json).map(([group, json]) => [group, JSON.parse(json)]));
  const stamped = Object.keys(options).length > 0;
  const raw = readFileSync(file), saved = kind === 'bin' ? await decodeState(new Uint8Array(raw)) : null;
  if (stamped && !saved) throw new Error('RADIATION and MOIST are recorded in binary states only');
  const ringed = !!saved && SST_FIELDS.some((name) => name in saved);
  const kept = stamped ? { ...withoutSstRing(saved), options: { ...saved.options, ...options } } : ringed ? withoutSstRing(saved) : null;
  const gz = kept ? gzipSync(encodeState(kept, { f64: f64Names(saved) }), { level: 9 }) : file.endsWith('.gz') ? raw : gzipSync(raw, { level: 9 });
  const parts = [];
  for (let at = 0; at < gz.length; at += chunk) {
    const name = `${base}.${kind}.gz.${String(parts.length).padStart(3, '0')}`;
    const bytes = gz.subarray(at, Math.min(at + chunk, gz.length));
    writeFileSync(name, bytes);
    parts.push({ file: name.replace(/.*\//, ''), bytes: bytes.length, sha256: sha256(bytes) });
  }
  writeFileSync(`${base}.parts.json`, manifest(base.replace(/.*\//, ''), null, parts));
  console.log(`${base}.parts.json: ${parts.length} parts, ${(gz.length / 1048576).toFixed(1)} MiB gzipped${ringed ? ', its SST ring left out' : ''}${stamped ? `, with the options ${JSON.stringify(kept.options)}` : ''}`);
}
