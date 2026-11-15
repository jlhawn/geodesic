import { LAKE_FIELDS, holdsLake, savedLake, savedLakeOptions } from '../js/physics/lake.module.js';

/*
 * The land of a saved state onto a model whose arrays are loaded; the deep
 * soil layers of a state saved without them start from the state's skin at
 * its time.
 */
export function loadSavedLand(model, saved) {
  model.land.load(saved.land, model.state[6], saved.time, model.state[3]);
}

export const savedIceOptions = (saved) => ({ ...(saved.given?.ocean?.iceExchange !== undefined ? { iceExchange: saved.given.ocean.iceExchange } : {}), ...(saved.given?.ice ?? {}) });

export function restoreSurface(model, saved, report = (line) => console.error(line)) {
  const { state, seaIce, land } = model;
  seaIce.load(state[6], saved.concentration ?? null);
  const exchange = seaIce.turbulent ? [['iceWater', seaIce.waterHeat], ['iceConductance', seaIce.conductance]].filter(([name, field]) => {
    const held = saved.ocean?.[name];
    if (held && held.length === field.length) field.set(held);
    return !(held && held.length === field.length);
  }).map(([name]) => `ocean.${name}`) : [];
  if (land && saved.land) land.load(saved.land, state[6], saved.time, state[3]);
  const missing = [...(saved.concentration ? [] : ['concentration']), ...exchange, ...(land && !saved.land ? ['land'] : []),
    ...(land && saved.land ? Object.keys(land.serialize()).filter((name) => !name.startsWith('goal') && !(name in saved.land)).map((name) => `land.${name}`) : [])];
  if (missing.length) report(`restoreSurface: the save holds no ${missing.join(', ')}; the model's starting values stand in for them`);
  return missing;
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
export function loadSavedLake(model, saved, { quiet = false } = {}) {
  const columns = savedLake(saved);
  if (!holdsLake(columns)) {
    if (model.lake?.cells.length && !quiet) console.error(`the save's lake columns are absent or never stepped: the ${model.lake.cells.length} lake cells start from fresh columns at their first step`);
    return model.lake && model.lake.read ? Promise.resolve(0) : 0;
  }
  if (!model.lake?.cells.length) throw new Error('the save holds lakes and the model has none');
  model.lake.load(columns);
  return model.lake.read ? model.lake.read().then(() => checkSavedLake(model.lake, columns)) : checkSavedLake(model.lake, columns);
}
