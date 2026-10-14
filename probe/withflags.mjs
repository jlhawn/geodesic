import { createRequire } from 'node:module';
const require = createRequire('/Users/jlhawn/git_repos/jlhawn/geodesic/package.json');
const { create, globals } = await import(require.resolve('webgpu'));
Object.assign(globalThis, globals);
const flags = (process.env.FLAGS ?? '').split(' ').filter(Boolean);
globalThis.navigator = { gpu: create(flags) };
console.log('dawn flags:', flags.join(' ') || '(none)');
await import('./bisect.mjs');
