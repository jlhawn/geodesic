import { readFileSync, realpathSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';

/*
 * The names a script reads from the environment are the ones its text
 * shows: every process.env.NAME and every quoted name in capitals, in the
 * script and in the helpers it imports from scripts/ that are not scripts
 * of their own (they read no process.argv). A name built at run time is
 * not seen, so an option group read under one is refused.
 */
export function envNames(url, names = new Set(), seen = new Set()) {
  const path = fileURLToPath(url);
  if (seen.has(path)) return names;
  seen.add(path);
  const text = readFileSync(path, 'utf8');
  for (const [, name] of text.matchAll(/process\.env\.([A-Za-z_$][\w$]*)/g)) names.add(name);
  for (const [, , name] of text.matchAll(/(['"`])([A-Z][A-Z0-9_]*)\1/g)) names.add(name);
  for (const [, from] of text.matchAll(/\bfrom\s*['"](\.{1,2}\/[^'"]+\.mjs)['"]/g)) {
    const helper = new URL(from, pathToFileURL(path));
    if (!/process\.argv/.test(readFileSync(helper, 'utf8'))) envNames(helper, names, seen);
  }
  return names;
}

export const unreadOptions = (url, env = process.env) => {
  const read = envNames(url);
  return Object.keys(env).filter((name) => /^\s*\{/.test(env[name]) && !read.has(name)).sort();
};

export function refuseUnreadOptions(url) {
  let main;
  try { main = realpathSync(process.argv[1]); } catch { return; }
  if (main !== fileURLToPath(url)) return;
  const unread = unreadOptions(url);
  if (unread.length) throw new Error(`${main.split('/').pop()} does not read ${unread.join(', ')}, given ${unread.length > 1 ? 'as option groups' : 'as an option group'}: unset or rename ${unread.length > 1 ? 'them' : 'it'}`);
}
