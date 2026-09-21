/*
 * A saved state is a JSON file, or that file gzipped and cut into parts
 * listed by a <name>.parts.json manifest so each piece stays under the
 * file-size caps of static hosts. Gzip is detected by its magic bytes
 * rather than the file name, since a server may already have inflated
 * a .gz file on the way.
 */
export function stateName(file) {
  return file.replace(/(?:\.json(?:\.gz)?|\.parts\.json)$/, '');
}

export function decodeState(bytes) {
  const gzipped = bytes.length > 1 && bytes[0] === 0x1f && bytes[1] === 0x8b;
  const body = new Response(bytes).body;
  return new Response(gzipped ? body.pipeThrough(new DecompressionStream('gzip')) : body).json();
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
