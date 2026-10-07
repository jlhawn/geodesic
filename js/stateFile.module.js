import { sigmaInterfaces, sigmaGridName } from './dynamics/sigmaCore.module.js';

/*
 * A saved state is a JSON file, or that file gzipped and cut into parts
 * listed by a <name>.parts.json manifest so each piece stays under the
 * file-size caps of static hosts, each part listed with its byte count
 * and, when the manifest was written with them, its SHA-256, the manifest
 * carrying the run's `id`, the SHA-256 of those hashes one per line (runId),
 * so that a run is named by its content and a store can keep its chunks
 * under their hashes; or the same state in binary: 'GCMS', a
 * little-endian u32 header length, the JSON header { N, K, day, time,
 * terrain, arrays: [{ name, length, offset, type }] } (which may end in
 * spaces counted in that length), zero padding to a multiple of 4, then
 * each array at its byte offset from the end of the padding, as float32
 * unless its type is 'f64'. Array names are the JSON
 * state's keys, with 'ocean.' and 'land.' prefixes for the nested ones.
 * Gzip and the binary format are recognised by their magic bytes rather
 * than the file name, since a server may already have inflated a .gz
 * file on the way. The array 'levels' holds the interface σ values of the
 * atmosphere's grid, top to ground; a state without it is on 'cam26'.
 */
export function stateName(file) {
  return file.replace(/(?:\.json(?:\.gz)?|\.parts\.json|\.bin(?:\.gz)?)$/, '');
}

/*
 * The saved states a directory index links to: the JSON, gzipped and
 * parts forms of <name>_state_day<N>, and the binary, gzipped and parts
 * forms of <name>_day<N>.
 */
export function listedStates(html) {
  return [...new Set([...html.matchAll(/href="([^"]+(?:_state_day\d+(?:\.json(?:\.gz)?|\.parts\.json)|_day\d+(?:\.bin(?:\.gz)?|\.parts\.json)))"/g)].map((m) => decodeURIComponent(m[1])))].sort();
}
export function stateDay(file) {
  const match = file.match(/_(?:state_)?day(\d+)/);
  return match ? Number(match[1]) : null;
}

const MAGIC = 'GCMS';
const isBinary = (bytes) => bytes.length > 8 && String.fromCharCode(bytes[0], bytes[1], bytes[2], bytes[3]) === MAGIC;

export async function decodeState(bytes) {
  const gzipped = bytes.length > 1 && bytes[0] === 0x1f && bytes[1] === 0x8b;
  if (!gzipped) return isBinary(bytes) ? decodeBinary(bytes) : new Response(bytes).json();
  const inflated = new Uint8Array(await new Response(new Response(bytes).body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
  return isBinary(inflated) ? decodeBinary(inflated) : JSON.parse(new TextDecoder().decode(inflated));
}

function decodeBinary(bytes) {
  const headerLength = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength).getUint32(4, true);
  const { arrays, ...state } = JSON.parse(new TextDecoder().decode(bytes.subarray(8, 8 + headerLength)));
  const start = 8 + Math.ceil(headerLength / 4) * 4;
  for (const { name, length, offset, type = 'f32' } of arrays) {
    const Type = type === 'f64' ? Float64Array : Float32Array, from = start + offset, to = from + length * Type.BYTES_PER_ELEMENT;
    if (to > bytes.byteLength) throw new Error(`${name} runs past the end of the state`);
    const values = (bytes.byteOffset + from) % Type.BYTES_PER_ELEMENT === 0 ? new Type(bytes.buffer, bytes.byteOffset + from, length) : new Type(bytes.slice(from, to).buffer);
    const [group, key] = name.includes('.') ? name.split('.') : [null, name];
    if (group) (state[group] ??= {})[key] = values; else state[key] = values;
  }
  return state;
}

/*
 * The sigma interfaces a saved state's atmosphere is on: those of the
 * named grid they match, else the saved values themselves.
 */
export function savedLevels(saved) {
  if (!saved || !saved.levels) return sigmaInterfaces('cam26');
  const name = sigmaGridName(saved.levels);
  return name ? sigmaInterfaces(name) : Float64Array.from(saved.levels);
}

/*
 * The binary form of a state: its top-level arrays and those under
 * `ocean` and `land`, as float32 or, listed in `f64` or named 'levels',
 * 'ocean.densities', 'land.record' or 'energyRecord', float64, with the header padded with spaces so the
 * arrays start 8-byte aligned.
 */
export function encodeState(state, { f64 = [] } = {}) {
  const named = [];
  for (const [key, value] of Object.entries(state)) {
    if (ArrayBuffer.isView(value) || Array.isArray(value)) named.push([key, value]);
    else if ((key === 'ocean' || key === 'land') && value) for (const [inner, values] of Object.entries(value)) named.push([`${key}.${inner}`, values]);
  }
  let offset = 0;
  const arrays = named.map(([name, values]) => {
    const type = f64.includes(name) || name === 'levels' || name === 'ocean.densities' || name === 'land.record' || name === 'energyRecord' ? 'f64' : 'f32', size = type === 'f64' ? 8 : 4;
    offset = Math.ceil(offset / size) * size;
    const entry = { name, length: values.length, offset, type };
    offset += size * values.length;
    return entry;
  });
  const meta = Object.fromEntries(Object.entries(state).filter(([, value]) => typeof value !== 'object'));
  const json = JSON.stringify({ ...meta, arrays });
  const header = new TextEncoder().encode(json.padEnd(json.length + ((8 - ((8 + new TextEncoder().encode(json).length) % 8)) % 8)));
  const start = 8 + header.length;
  const bytes = new Uint8Array(start + offset);
  for (let c = 0; c < 4; c++) bytes[c] = MAGIC.charCodeAt(c);
  new DataView(bytes.buffer).setUint32(4, header.length, true);
  bytes.set(header, 8);
  arrays.forEach(({ offset: at, length, type }, n) => new (type === 'f64' ? Float64Array : Float32Array)(bytes.buffer, start + at, length).set(named[n][1]));
  return bytes;
}

async function sha256Hex(bytes) {
  const digest = await crypto.subtle.digest('SHA-256', bytes);
  return Array.from(new Uint8Array(digest), (b) => b.toString(16).padStart(2, '0')).join('');
}

// A run's id: the SHA-256 of its parts' hashes, one per line.
export function runId(parts) {
  return sha256Hex(new TextEncoder().encode(parts.map((part) => part.sha256).join('\n')));
}

/*
 * The parts a state's address names, with the run they make when every
 * part carries its hash. A manifest that cannot be fetched is taken from
 * the store, which has it once the run has loaded from there.
 */
async function manifestOf(url, store) {
  if (!/\.parts\.json(?:[?#]|$)/.test(url)) return { parts: [{ url, bytes: 0, sha256: null }], run: null };
  let manifest;
  try {
    const response = await fetch(url);
    if (!response.ok) throw new Error(`${url}: ${response.status}`);
    manifest = await response.json();
  } catch (error) {
    const kept = store ? await store.runForUrl(url).catch(() => null) : null;
    if (!kept) throw error;
    manifest = kept;
  }
  const parts = manifest.parts.map((part) => ({ url: new URL(part.file, url).href, file: part.file, bytes: part.bytes, sha256: part.sha256 ?? null }));
  if (!parts.length || !parts.every((part) => part.sha256)) return { parts, run: null };
  const name = manifest.name ?? stateName(url.replace(/.*\//, '').replace(/[?#].*/, ''));
  const run = { id: await runId(parts), name, description: manifest.description ?? null, url, parts: parts.map(({ file, bytes, sha256 }) => ({ file, bytes, sha256 })), bytes: parts.reduce((sum, part) => sum + part.bytes, 0) };
  return { parts, run };
}

function joined(pieces, length) {
  if (pieces.length === 1) return pieces[0];
  const bytes = new Uint8Array(length);
  let at = 0;
  for (const piece of pieces) { bytes.set(piece, at); at += piece.length; }
  return bytes;
}

async function* chunksOf(reader) {
  try {
    for (;;) {
      const { done, value } = await reader.read();
      if (done) return;
      yield value;
    }
  } finally {
    reader.cancel().catch(() => {});
  }
}

async function* inflated(chunks) {
  const source = new ReadableStream({
    async pull(controller) {
      const { done, value } = await chunks.next();
      if (done) controller.close(); else controller.enqueue(value);
    },
    async cancel() { await chunks.return(); },
  });
  yield* chunksOf(source.pipeThrough(new DecompressionStream('gzip')).getReader());
}

async function* following(head, rest) {
  if (head.length) yield head;
  yield* rest;
}

async function gather(chunks, wanted, head = new Uint8Array(0)) {
  const pieces = [head];
  let length = head.length;
  while (length < wanted) {
    const { done, value } = await chunks.next();
    if (done) break;
    pieces.push(value);
    length += value.length;
  }
  if (pieces.length === 1) return head;
  if (pieces.length === 2 && !head.length) return pieces[1];
  const bytes = new Uint8Array(length);
  let at = 0;
  for (const piece of pieces) { bytes.set(piece, at); at += piece.length; }
  return bytes;
}

/*
 * A state's parts are fetched in turn as one stream of bytes, inflated on
 * the way when they are gzipped. A binary state goes into one buffer of
 * the size its header gives, which decodeBinary reads in place; any other
 * state is gathered and decoded as decodeState decodes it. With a store
 * (runStore in snapshots.module.js, or anything with its methods), a part
 * that names its hash comes from the store when it is there and otherwise
 * from the network, checked against the hash and put in the store, and a
 * run whose every part did so is recorded there once it has decoded, so
 * the next load has it without the network; a store that fails only
 * leaves the parts on the network. `progress` gets the bytes so far, the
 * total and whether the part came from the 'store' or the 'network'.
 */
export async function fetchState(url, progress = null, { store = null } = {}) {
  const { parts, run } = await manifestOf(url, store);
  let total = parts.reduce((sum, part) => sum + part.bytes, 0), received = 0, keeping = !!(store && run);
  async function* fetched() {
    for (const part of parts) {
      const kept = store && part.sha256 ? await store.getChunk(part.sha256).catch(() => null) : null;
      if (kept) {
        received += kept.length;
        progress?.(received, total, 'store');
        yield kept;
        continue;
      }
      const response = await fetch(part.url);
      if (!response.ok) throw new Error(`${part.url}: ${response.status}`);
      if (!total) total = Number(response.headers.get('content-length')) || 0;
      const body = response.body ? chunksOf(response.body.getReader()) : [new Uint8Array(await response.arrayBuffer())];
      if (!part.sha256) {
        for await (const value of body) {
          received += value.length;
          progress?.(received, total, 'network');
          yield value;
        }
        continue;
      }
      const pieces = [];
      let length = 0;
      for await (const value of body) {
        pieces.push(value);
        length += value.length;
        received += value.length;
        progress?.(received, total, 'network');
      }
      const bytes = joined(pieces, length), digest = await sha256Hex(bytes);
      if (digest !== part.sha256) throw new Error(`${part.file}: the part's SHA-256 is ${digest.slice(0, 12)}…, not the ${part.sha256.slice(0, 12)}… its manifest names`);
      if (keeping) keeping = await store.putChunk(part.sha256, bytes).then(() => true, () => false);
      yield bytes;
    }
  }
  const state = await decoded(url, fetched());
  if (keeping) await store.putRun(run).then(() => store.trimRuns()).catch(() => {});
  return state;
}

async function decoded(url, raw) {
  const first = await gather(raw, 2);
  const chunks = first.length > 1 && first[0] === 0x1f && first[1] === 0x8b ? inflated(following(first, raw)) : following(first, raw);
  let head = await gather(chunks, 9);
  if (!isBinary(head)) return decodeState(await gather(chunks, Infinity, head));
  const headerLength = new DataView(head.buffer, head.byteOffset, head.byteLength).getUint32(4, true);
  head = await gather(chunks, 8 + headerLength, head);
  const { arrays } = JSON.parse(new TextDecoder().decode(head.subarray(8, 8 + headerLength)));
  const start = 8 + Math.ceil(headerLength / 4) * 4;
  const ends = arrays.map(({ name, length, offset, type = 'f32' }) => ({ name, end: start + offset + length * (type === 'f64' ? 8 : 4) }));
  const bytes = new Uint8Array(ends.reduce((size, { end }) => Math.max(size, end), start));
  if (head.length > bytes.length) throw new Error(`${url}: the state runs past the end of its arrays`);
  bytes.set(head);
  let at = head.length;
  for await (const chunk of chunks) {
    if (at + chunk.length > bytes.length) throw new Error(`${url}: the state runs past the end of its arrays`);
    bytes.set(chunk, at);
    at += chunk.length;
  }
  const short = ends.find(({ end }) => end > at);
  if (short) throw new Error(`${short.name} runs past the end of the state`);
  return decodeBinary(bytes);
}
