/*
 * When the address names no resolution, engine or saved run, the worker's
 * device test (probe in model.worker.js) picks them: the highest
 * resolution whose projected rate clears TARGET_RATE simulated hours a
 * minute, 24 plus room for the page's own work. A step's time grows with
 * the number of cells, N², from the GPU's test at N=64, or from one CPU
 * thread's at N=16 divided by the workers' speed-up. The choice is kept in
 * this browser for the same browser and GPU.
 */
export const TARGET_RATE = 30, PROBE_VERSION = 1, GPU_LADDER = [128, 64, 32], CPU_LADDER = [64, 32, 16];
export const hoursPerMinute = (N, ms) => (1350 * 16 / N) / (ms / 1000) / 60;
export function pickDevice({ gpu, cpu }, workers) {
  const pick = (ladder, test, speedup) => {
    const rates = ladder.map((N) => hoursPerMinute(N, test.ms * (N / test.N) ** 2 / speedup));
    const first = rates.findIndex((rate) => rate >= TARGET_RATE), at = first < 0 ? ladder.length - 1 : first;
    return { N: ladder[at], rate: rates[at] };
  };
  if (gpu && !gpu.error) return { engine: 'gpu', ...pick(GPU_LADDER, gpu, 1), measured: `${gpu.ms.toFixed(1)} ms a step at N=${gpu.N} on the GPU` };
  if (!cpu || cpu.error) return null;
  const speedup = workers > 1 ? Math.min(0.7 * workers, 7) : 1;
  return { engine: 'cpu', ...pick(CPU_LADDER, cpu, speedup), measured: `${cpu.ms.toFixed(0)} ms a step at N=${cpu.N} on one CPU thread${gpu ? `; the GPU failed: ${gpu.error}` : ''}` };
}
