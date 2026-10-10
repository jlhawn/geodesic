export async function loadFrames(url) {
  const bytes = await (await fetch(url)).arrayBuffer();
  const head = new DataView(bytes);
  if (String.fromCharCode(...new Uint8Array(bytes, 0, 4)) !== 'EXF1') throw new Error(`${url} is not a frames file`);
  const length = head.getUint32(4, true), header = JSON.parse(new TextDecoder().decode(new Uint8Array(bytes, 8, length)));
  const C = 10 * header.N * header.N + 2, F = header.fields.length;
  let offset = 8 + length;
  const sequences = {};
  for (const s of header.sequences) { sequences[s.name] = { times: s.times, start: offset }; offset += s.times.length * F * C * 2; }
  return {
    header, C,
    values(sequence, frame, field, out = new Float32Array(C)) {
      const { scale, offset: base } = header.fields[field], raw = new Int16Array(bytes, sequences[sequence].start + (frame * F + field) * C * 2, C);
      for (let i = 0; i < C; i++) out[i] = raw[i] * scale + base;
      return out;
    },
    times: (sequence) => sequences[sequence].times,
  };
}
