// The files of a run's OUT directory that the Mac pulls from a Verda
// instance, listed the same way on both sides so that scripts/verdaPull.sh
// can compare them: one line per file, sorted by path,
//   <size in bytes>\t<sha256 hex, or '-' without --sha256>\t<path under OUT>
// Every regular file under OUT counts except files still being written
// (*.partial), in-day checkpoints (*_stepSSSS.bin), lock directories and
// any other .bin than a whole-day state of PREFIX at the top level
// (<PREFIX><N>_dayDDDD.bin), so the benchmark's states stay behind.
//   node scripts/verdaFiles.mjs <OUT> <PREFIX> [--sha256]
import { createReadStream, readdirSync, statSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

export function wanted(path, prefix) {
  const parts = path.split('/'), name = parts[parts.length - 1];
  if (parts.some((part) => part.endsWith('.lock'))) return false;
  if (name.endsWith('.partial')) return false;
  if (name.endsWith('.bin')) return parts.length === 1 && new RegExp(`^${prefix}\\d+_day\\d+\\.bin$`).test(name);
  return true;
}

export function listFiles(dir, prefix, at = '') {
  const found = [];
  for (const entry of readdirSync(join(dir, at), { withFileTypes: true })) {
    const path = at ? `${at}/${entry.name}` : entry.name;
    if (entry.isDirectory()) { if (!entry.name.endsWith('.lock')) found.push(...listFiles(dir, prefix, path)); }
    else if (entry.isFile() && wanted(path, prefix)) found.push(path);
  }
  return found.sort();
}

const sha256 = (file) => new Promise((resolve, reject) => {
  const hash = createHash('sha256');
  createReadStream(file).on('data', (chunk) => hash.update(chunk)).on('end', () => resolve(hash.digest('hex'))).on('error', reject);
});

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const args = process.argv.slice(2), withSums = args.includes('--sha256');
  const [dir, prefix] = args.filter((a) => a !== '--sha256');
  if (!dir || !prefix) { console.error('usage: node scripts/verdaFiles.mjs <OUT> <PREFIX> [--sha256]'); process.exit(1); }
  for (const path of listFiles(dir, prefix)) {
    const file = join(dir, path);
    console.log(`${statSync(file).size}\t${withSums ? await sha256(file) : '-'}\t${path}`);
  }
}
