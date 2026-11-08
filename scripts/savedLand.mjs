import { LAKE_FIELDS, holdsLake, savedLake, savedLakeOptions } from '../js/physics/lake.module.js';

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

export const lakeOptionsFor = (saved) => ({ ...savedLakeOptions(saved), ...JSON.parse(process.env.LAKE ?? '{}') });

export function savedLakeOff(field, columns, cells) {
  return LAKE_FIELDS.filter((name) => columns[name]).map((name) => {
    const values = field(name), want = columns[name];
    return [name, cells.filter((i) => !(Math.abs(values[i] - want[i]) <= 1e-6 * Math.max(1, Math.abs(want[i])))).length];
  }).filter(([, count]) => count > 0);
}

function checkSavedLake(lake, columns) {
  const off = savedLakeOff((name) => lake[name], columns, lake.cells);
  if (off.length) throw new Error(`the lake columns differ from the save after the load (cells off): ${off.map(([name, count]) => `${name} ${count}`).join(', ')}`);
  return lake.cells.length;
}

/*
 * A saved state's lake columns onto a model at its resolution, checked
 * against the save: on the GPU model after reading them back from the
 * device, so the result is a promise there. Without them each lake starts
 * a fresh column at its first step and drops the ice the save holds.
 */
export function loadSavedLake(model, saved) {
  const columns = savedLake(saved);
  if (!holdsLake(columns)) {
    if (model.lake?.cells.length) console.error(`the save holds no lake columns: the ${model.lake.cells.length} lake cells start from fresh columns at their first step`);
    return model.lake && model.lake.read ? Promise.resolve(0) : 0;
  }
  if (!model.lake?.cells.length) throw new Error('the save holds lakes and the model has none');
  model.lake.load(columns);
  return model.lake.read ? model.lake.read().then(() => checkSavedLake(model.lake, columns)) : checkSavedLake(model.lake, columns);
}
