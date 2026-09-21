// Usage: node scripts/splitState.mjs runs/<name>_state_dayNNN.json [chunkMiB]
// Writes <name>.json.gz.000, .001, … of at most chunkMiB each and the
// <name>.parts.json manifest the page loads them from.
import { readFileSync, writeFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';

const [file, mib = '10'] = process.argv.slice(2);
if (!file) { console.error('usage: node scripts/splitState.mjs <state.json|state.json.gz> [chunkMiB]'); process.exit(1); }
const chunk = Number(mib) * 1048576;
const base = file.replace(/\.json(?:\.gz)?$/, '');
const gz = file.endsWith('.gz') ? readFileSync(file) : gzipSync(readFileSync(file), { level: 9 });
const parts = [];
for (let at = 0; at < gz.length; at += chunk) {
  const name = `${base}.json.gz.${String(parts.length).padStart(3, '0')}`;
  const bytes = gz.subarray(at, Math.min(at + chunk, gz.length));
  writeFileSync(name, bytes);
  parts.push({ file: name.replace(/.*\//, ''), bytes: bytes.length });
}
writeFileSync(`${base}.parts.json`, JSON.stringify({ name: base.replace(/.*\//, ''), parts }, null, 1) + '\n');
console.log(`${base}.parts.json: ${parts.length} parts, ${(gz.length / 1048576).toFixed(1)} MiB gzipped`);
