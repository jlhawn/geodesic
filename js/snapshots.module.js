/*
 * Browser-side snapshot store for the climate page, in IndexedDB: a
 * `meta` store with the name, creation time and size of each snapshot
 * and a `data` store with the state arrays, so listing never reads the
 * tens of megabytes of state.
 */
const DB_NAME = 'climate', VERSION = 1;

function openStore() {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(DB_NAME, VERSION);
    request.onupgradeneeded = () => {
      const db = request.result;
      if (!db.objectStoreNames.contains('meta')) db.createObjectStore('meta', { keyPath: 'id', autoIncrement: true });
      if (!db.objectStoreNames.contains('data')) db.createObjectStore('data', { keyPath: 'id' });
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
