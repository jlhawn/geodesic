/*
 * A saved state is a JSON file, or that file gzipped and cut into parts
 * listed by a <name>.parts.json manifest so each piece stays under the
 * file-size caps of static hosts, or the same state in binary: 'GCMS', a
 * little-endian u32 header length, the JSON header { N, K, day, time,
 * terrain, arrays: [{ name, length, offset, type }] }, zero padding to a
 * multiple of 4, then each array at its byte offset from the end of the
 * padding, as float32 unless its type is 'f64'. Array names are the JSON
 * state's keys, with 'ocean.' and 'land.' prefixes for the nested ones.
 * Gzip and the binary format are recognised by their magic bytes rather
 * than the file name, since a server may already have inflated a .gz
 * file on the way.
 */
export function stateName(file) {
  return file.replace(/(?:\.json(?:\.gz)?|\.parts\.json|\.bin(?:\.gz)?)$/, '');
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
 * The binary form of a state: its top-level arrays and those under
 * `ocean` and `land`, as float32 or, listed in `f64`, float64.
 */
export function encodeState(state, { f64 = [] } = {}) {
  const named = [];
  for (const [key, value] of Object.entries(state)) {
    if (ArrayBuffer.isView(value) || Array.isArray(value)) named.push([key, value]);
    else if ((key === 'ocean' || key === 'land') && value) for (const [inner, values] of Object.entries(value)) named.push([`${key}.${inner}`, values]);
  }
  let offset = 0;
  const arrays = named.map(([name, values]) => {
    const type = f64.includes(name) ? 'f64' : 'f32', size = type === 'f64' ? 8 : 4;
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

async function partsOf(url) {
  if (!/\.parts\.json(?:[?#]|$)/.test(url)) return [{ url, bytes: 0 }];
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: ${response.status}`);
  const manifest = await response.json();
  return manifest.parts.map((part) => ({ url: new URL(part.file, url).href, bytes: part.bytes }));
}

export async function fetchState(url, progress = null) {
  const parts = await partsOf(url);
  let total = parts.reduce((sum, part) => sum + part.bytes, 0);
  const chunks = [];
  let received = 0;
  for (const part of parts) {
    const response = await fetch(part.url);
    if (!response.ok) throw new Error(`${part.url}: ${response.status}`);
    if (!total) total = Number(response.headers.get('content-length')) || 0;
    if (!response.body) {
      const bytes = new Uint8Array(await response.arrayBuffer());
      chunks.push(bytes);
      received += bytes.length;
      progress?.(received, total);
      continue;
    }
    const reader = response.body.getReader();
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value);
      received += value.length;
      progress?.(received, total);
    }
  }
  const bytes = new Uint8Array(received);
  let at = 0;
  for (const chunk of chunks) { bytes.set(chunk, at); at += chunk.length; }
  return decodeState(bytes);
}
