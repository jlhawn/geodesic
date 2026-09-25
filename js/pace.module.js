/*
 * The idle pause the model worker takes between GPU steps, chosen from
 * the page's once-a-second count of frames that came late for want of
 * the GPU. A pause costs simulation speed, so it is worth taking only when
 * the model's GPU work is what makes frames late. After `raise` late
 * seconds in a row the pacer tries the next level for `trial` reports,
 * and keeps it only if the late frames fell to at most `keep` of their
 * rate at the level below. A level that didn't help is not tried again
 * for `cooldown` reports. After `calm` quiet reports (fewer than
 * `threshold` late frames) the pacer steps back down a level. A step-down
 * that lateness reverses at once doubles the calm needed before the next,
 * up to `calmMax`, and `calmMax` quiet reports without a pause restore it.
 */
export const PAUSE_LEVELS = [0, 2, 4, 8];

export function createPacer({ raise = 2, trial = 3, keep = 2 / 3, cooldown = 60, calm = 5, calmMax = 80, threshold = 2 } = {}) {
  let level = 0, reports = 0, lateRun = 0, quiet = 0, calmNeeded = calm, downAt = -Infinity;
  let trying = null, blockedFrom = Infinity, blockedUntil = 0;
  const recent = [];
  const mean = (values) => values.reduce((sum, v) => sum + v, 0) / values.length;

  function report(late) {
    reports++;
    recent.push(late);
    if (recent.length > trial) recent.shift();
    if (trying) {
      trying.seen.push(late);
      if (trying.seen.length < trial) return PAUSE_LEVELS[level];
      if (mean(trying.seen) <= keep * mean(trying.before)) {
        if (reports - downAt <= raise + trial + 1) calmNeeded = Math.min(calmMax, 2 * calmNeeded);
      } else {
        level = trying.from;
        blockedFrom = trying.from + 1;
        blockedUntil = reports + cooldown;
      }
      trying = null;
      recent.length = 0;
      return PAUSE_LEVELS[level];
    }
    if (late >= threshold) {
      quiet = 0;
      const blocked = level + 1 >= blockedFrom && reports < blockedUntil;
      if (++lateRun >= raise && level + 1 < PAUSE_LEVELS.length && !blocked) {
        trying = { from: level, before: recent.slice(), seen: [] };
        level++;
        lateRun = 0;
      }
    } else {
      lateRun = 0;
      quiet++;
      if (level > 0 && quiet >= calmNeeded) { level--; quiet = 0; downAt = reports; }
      else if (level === 0 && quiet >= calmMax) calmNeeded = calm;
    }
    return PAUSE_LEVELS[level];
  }

  return { report, pause: () => PAUSE_LEVELS[level] };
}
