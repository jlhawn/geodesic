/*
 * The thread primitives the parallel engine needs, from worker_threads in
 * Node and from Web Workers in the browser: spawning a worker with an
 * initial message, and receiving that message inside the worker.
 */
const inBrowser = typeof process === 'undefined' || !process.versions || !process.versions.node;

export async function parallelism() {
  if (inBrowser) return navigator.hardwareConcurrency || 4;
  return (await import('node:os')).availableParallelism();
}

export async function spawn(url, data, { onMessage, onError }) {
  if (inBrowser) {
    const worker = new Worker(url, { type: 'module' });
    worker.onmessage = (event) => onMessage(event.data);
    worker.onerror = (event) => onError(event.error ?? new Error(event.message || `${url} failed to load (${event.filename}:${event.lineno})`));
    worker.postMessage(data);
    return { terminate: () => worker.terminate(), unref() {} };
  }
  const { Worker: Thread } = await import('node:worker_threads');
  const thread = new Thread(url, { workerData: data });
  thread.on('message', onMessage);
  thread.on('error', onError);
  return { terminate: () => thread.terminate(), unref: () => thread.unref() };
}

export async function workerInit() {
  if (inBrowser) {
    const data = await new Promise((resolve) => { self.onmessage = (event) => resolve(event.data); });
    return { data, post: (message) => self.postMessage(message) };
  }
  const { parentPort, workerData } = await import('node:worker_threads');
  return { data: workerData, post: (message) => parentPort.postMessage(message) };
}
