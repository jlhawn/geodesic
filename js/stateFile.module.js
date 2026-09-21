/*
 * A saved state is JSON, gzipped when it has to fit GitHub's file
 * limits; the gzip magic bytes decide rather than the file name, since a
 * server may already have inflated a .gz file on the way.
 */
export function stateName(file) {
  return file.replace(/\.json(?:\.gz)?$/, '');
}

export function decodeState(bytes) {
  const gzipped = bytes.length > 1 && bytes[0] === 0x1f && bytes[1] === 0x8b;
  const body = new Response(bytes).body;
  return new Response(gzipped ? body.pipeThrough(new DecompressionStream('gzip')) : body).json();
}

export async function readState(response) {
  if (!response.ok) throw new Error(`${response.url}: ${response.status}`);
  return decodeState(new Uint8Array(await response.arrayBuffer()));
}
