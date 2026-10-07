/*
 * Browser-side stores for the climate page, in IndexedDB. The custom
 * snapshots: a `meta` store with the name, creation time and size of each
 * and a `data` store with the state arrays, so listing never reads the
 * tens of megabytes of state. The saved runs the page has fetched: their
 * compressed chunks in `chunks` under their SHA-256 and each run's
 * manifest in `runs` under the run's id, so a run loads from here once it
 * has been downloaded, runs that share chunks share the bytes, and the
 * least recently used runs go once the chunks pass RUN_STORE_LIMIT.
 */
const DB_NAME = 'climate', VERSION = 2;
export const RUN_STORE_LIMIT = 1024 * 1048576;
const ORPHAN_AGE = 3600e3;

function openStore() {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(DB_NAME, VERSION);
    request.onupgradeneeded = () => {
      const db = request.result;
      if (!db.objectStoreNames.contains('meta')) db.createObjectStore('meta', { keyPath: 'id', autoIncrement: true });
      if (!db.objectStoreNames.contains('data')) db.createObjectStore('data', { keyPath: 'id' });
      if (!db.objectStoreNames.contains('chunks')) db.createObjectStore('chunks', { keyPath: 'sha256' }).createIndex('stored', 'stored');
      if (!db.objectStoreNames.contains('runs')) db.createObjectStore('runs', { keyPath: 'id' }).createIndex('url', 'url');
    };
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
}

function transaction(db, stores, mode, work) {
  return new Promise((resolve, reject) => {
    const tx = db.transaction(stores, mode);
    let result;
    tx.oncomplete = () => resolve(result);
    tx.onerror = () => reject(tx.error);
    tx.onabort = () => reject(tx.error);
    work(tx, (value) => { result = value; });
  });
}

const await1 = (request) => new Promise((resolve, reject) => { request.onsuccess = () => resolve(request.result); request.onerror = () => reject(request.error); });

export async function listSnapshots() {
  const db = await openStore();
  const all = await transaction(db, ['meta'], 'readonly', (tx, done) => { tx.objectStore('meta').getAll().onsuccess = (e) => done(e.target.result); });
  db.close();
  return all.sort((a, b) => b.created - a.created);
}

export async function saveSnapshot(meta, data) {
  const db = await openStore();
  const id = await transaction(db, ['meta', 'data'], 'readwrite', (tx, done) => {
    tx.objectStore('meta').add(meta).onsuccess = (e) => { const id = e.target.result; tx.objectStore('data').put({ id, ...data }); done(id); };
  });
  db.close();
  return id;
}

export async function getSnapshot(id) {
  const db = await openStore();
  const tx = db.transaction(['meta', 'data'], 'readonly');
  const [meta, data] = await Promise.all([await1(tx.objectStore('meta').get(id)), await1(tx.objectStore('data').get(id))]);
  db.close();
  return { meta, data };
}

export async function renameSnapshot(id, name) {
  const db = await openStore();
  await transaction(db, ['meta'], 'readwrite', (tx) => {
    const store = tx.objectStore('meta');
    store.get(id).onsuccess = (e) => { const meta = e.target.result; if (meta) store.put({ ...meta, name }); };
  });
  db.close();
}

export async function deleteSnapshot(id) {
  const db = await openStore();
  await transaction(db, ['meta', 'data'], 'readwrite', (tx) => { tx.objectStore('meta').delete(id); tx.objectStore('data').delete(id); });
  db.close();
}

export async function cloneSnapshot(id, name) {
  const { meta, data } = await getSnapshot(id);
  const copy = { ...meta, name, created: Date.now() };
  delete copy.id;
  const { id: _, ...payload } = data;
  return saveSnapshot(copy, payload);
}

export async function getChunk(sha256) {
  const db = await openStore();
  const row = await transaction(db, ['chunks'], 'readonly', (tx, done) => { tx.objectStore('chunks').get(sha256).onsuccess = (e) => done(e.target.result); });
  db.close();
  return row ? new Uint8Array(row.bytes) : null;
}

export async function putChunk(sha256, bytes) {
  const whole = bytes.byteOffset === 0 && bytes.byteLength === bytes.buffer.byteLength ? bytes.buffer : bytes.slice().buffer;
  const db = await openStore();
  await transaction(db, ['chunks'], 'readwrite', (tx) => { tx.objectStore('chunks').put({ sha256, bytes: whole, stored: Date.now() }); });
  db.close();
}

export async function listRuns() {
  const db = await openStore();
  const all = await transaction(db, ['runs'], 'readonly', (tx, done) => { tx.objectStore('runs').getAll().onsuccess = (e) => done(e.target.result); });
  db.close();
  return all.sort((a, b) => b.used - a.used);
}

export async function runForUrl(url) {
  const db = await openStore();
  const all = await transaction(db, ['runs'], 'readonly', (tx, done) => { tx.objectStore('runs').index('url').getAll(url).onsuccess = (e) => done(e.target.result); });
  db.close();
  return all.sort((a, b) => b.used - a.used)[0] ?? null;
}

export async function putRun(run) {
  const db = await openStore();
  await transaction(db, ['runs'], 'readwrite', (tx) => {
    const store = tx.objectStore('runs');
    store.get(run.id).onsuccess = (e) => { const now = Date.now(); store.put({ ...run, stored: e.target.result?.stored ?? now, used: now }); };
  });
  db.close();
}

const hashesOf = (runs) => new Set(runs.flatMap((run) => run.parts.map((part) => part.sha256)));
export const runBytes = (runs) => { const seen = new Map(); for (const run of runs) for (const part of run.parts) seen.set(part.sha256, part.bytes); return [...seen.values()].reduce((sum, bytes) => sum + bytes, 0); };

function dropRun(tx, runs, gone, kept) {
  const chunks = tx.objectStore('chunks'), shared = hashesOf(kept);
  for (const part of gone.parts) if (!shared.has(part.sha256)) chunks.delete(part.sha256);
  runs.delete(gone.id);
}

export async function deleteRun(id) {
  const db = await openStore();
  await transaction(db, ['runs', 'chunks'], 'readwrite', (tx) => {
    const runs = tx.objectStore('runs');
    runs.getAll().onsuccess = (e) => {
      const all = e.target.result, gone = all.find((run) => run.id === id);
      if (gone) dropRun(tx, runs, gone, all.filter((run) => run !== gone));
    };
  });
  db.close();
}

/*
 * Drops the least recently used runs while the chunks they share pass the
 * limit, the most recent always staying, then any chunk no run names that
 * has sat for ORPHAN_AGE, a download that never finished.
 */
export async function trimRuns(limit = RUN_STORE_LIMIT) {
  const db = await openStore();
  await transaction(db, ['runs', 'chunks'], 'readwrite', (tx) => {
    const runs = tx.objectStore('runs'), chunks = tx.objectStore('chunks');
    runs.getAll().onsuccess = (e) => {
      const all = e.target.result.sort((a, b) => b.used - a.used);
      while (all.length > 1 && runBytes(all) > limit) dropRun(tx, runs, all.pop(), all);
      const named = hashesOf(all);
      chunks.index('stored').openKeyCursor(IDBKeyRange.upperBound(Date.now() - ORPHAN_AGE)).onsuccess = (ev) => {
        const cursor = ev.target.result;
        if (!cursor) return;
        if (!named.has(cursor.primaryKey)) chunks.delete(cursor.primaryKey);
        cursor.continue();
      };
    };
  });
  db.close();
}

export const runStore = { getChunk, putChunk, runForUrl, putRun, trimRuns };
