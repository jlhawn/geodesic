// The warm layer at local 03 h, counted only where a warm layer there is not physical: in cells that are dark from local 00 h to 03 h and
// whose net surface heat flux over those hours was negative. Under polar day the sun is up, and under warm air over cold water the sea
// gains heat; in both a warm layer at 03 h is physical.
//   const census = createNightCensus(C);
//   census.observe(i, hour, insolation, netFlux, dt, warm)   for every open-water cell after each physics pass (hour is local solar time)
//   census.take()  { dark, darkWarm, darkOver01, darkOver1, darkWarmest, darkCell } since the last take, then resets the counts
// A cell is sampled in local hours 2.5 to 3.5. The window since 00 h restarts when a cell's hour wraps from 24 to 0.
export function createNightCensus(C) {
  const last = new Float64Array(C).fill(24), lit = new Uint8Array(C), net = new Float64Array(C);
  let t = null;
  const reset = () => { t = { dark: 0, darkWarm: 0, darkOver01: 0, darkOver1: 0, darkWarmest: 0, darkCell: -1 }; };
  reset();
  return {
    observe(i, hour, insolation, netFlux, dt, warm) {
      if (hour < last[i]) { lit[i] = 0; net[i] = 0; }
      last[i] = hour;
      if (hour >= 3.5) return;
      if (insolation > 0) lit[i] = 1;
      net[i] += netFlux * dt;
      if (hour < 2.5 || lit[i] || !(net[i] < 0)) return;
      t.dark++;
      if (!(warm > 0)) return;
      t.darkWarm++;
      if (warm > 0.1) t.darkOver01++;
      if (warm > 1) t.darkOver1++;
      if (warm > t.darkWarmest) { t.darkWarmest = warm; t.darkCell = i; }
    },
    take() { const out = t; reset(); return out; },
  };
}
