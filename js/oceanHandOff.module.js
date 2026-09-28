/*
 * The hand-off of asynchronous coupling (scripts/asyncSpinup.sh): a
 * decoded coupled state whose ocean, sea ice and sea surface are taken
 * from another decoded state of the same mesh, typically the end of an
 * ocean-only spin-up (scripts/oceanSpinup.mjs) started from an earlier
 * state of the same coupled run. Over the sea (where `land` is 0) the
 * surface temperature, ice, concentration and the snow on the ice come
 * from `ocean`, as do all the ocean's layers; the atmosphere with the
 * mixed-layer deck's carried state (mlmSubsidence, mlmHeight, mlmGate),
 * the land cells, the day and the time stay those of `coupled`.
 * `oceanYears`, the model years the ocean has spent alone, adds up the
 * two.
 */
export function withOceanOf(coupled, ocean, land) {
  if (coupled.N !== ocean.N) throw new Error(`the ocean comes from N=${ocean.N}, the coupled state is N=${coupled.N}`);
  const C = coupled.surfaceT.length;
  if (land.length !== C || ocean.surfaceT.length !== C) throw new Error('the land mask and both states must have one value a cell');
  if (!ocean.ocean || !ocean.ocean.h) throw new Error('the state to take the ocean from has no ocean');
  const sea = (mine, theirs) => Float32Array.from(mine, (x, i) => (land[i] ? x : theirs[i]));
  const concentration = ocean.concentration ? (coupled.concentration ? sea(coupled.concentration, ocean.concentration) : Float32Array.from(ocean.concentration)) : null;
  return {
    ...coupled,
    surfaceT: sea(coupled.surfaceT, ocean.surfaceT),
    ice: sea(coupled.ice, ocean.ice),
    concentration,
    ocean: { ...ocean.ocean },
    land: { ...coupled.land, snow: sea(coupled.land.snow, ocean.land.snow) },
    oceanYears: (coupled.oceanYears ?? 0) + (ocean.oceanYears ?? 0),
  };
}
