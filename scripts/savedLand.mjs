/*
 * The land of a saved state onto a model whose arrays are loaded. The ice
 * argument is left to the model's own default, which on the GPU model is
 * its sea ice, so the snow a save carries on sea ice stays; the deep soil
 * layers of a state saved without them start from the state's skin at its
 * time.
 */
export function loadSavedLand(model, saved) {
  model.land.load(saved.land, undefined, saved.time, model.state[3]);
}
