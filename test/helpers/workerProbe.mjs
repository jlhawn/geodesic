// The page's worker answering a device test as it does in a browser: no node:fs in the modules, every file through fetch.
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

globalThis.fetch = async (input) => {
  const path = fileURLToPath(new URL(typeof input === 'string' ? input : input.href ?? input.url));
  if (!existsSync(path)) return new Response(null, { status: 404 });
  return new Response(readFileSync(path), { status: 200 });
};
Object.defineProperty(process, 'getBuiltinModule', { value: undefined, configurable: true });
globalThis.self = { postMessage: (message) => { if (message.type === 'probe') { console.log(JSON.stringify(message)); process.exit(0); } } };
await import('../../js/model.worker.js');
await self.onmessage({ data: { type: 'probe', engine: 'cpu', land: true, terrain: true, topography: null } });
