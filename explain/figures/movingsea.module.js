import { playbackGlobe } from '../playbackGlobe.module.js';
import { loadFrames } from '../frames.module.js';
import { rampRGB, clamp, DARK_NEUTRAL, VAPOR } from '../runtime.module.js';
import { DATA_URL, calendarDate, season } from '../oceanRun.module.js';
import { seasonPhrase } from '../insolation.module.js';

const DARK = 'rgb(52, 55, 62)', PALE = '#d4efff', CURRENT = 'rgb(80, 200, 185)';
const hex = (color) => (color.startsWith('#') ? color.slice(1).match(/[0-9a-f]{2}/gi).map((h) => parseInt(h, 16)) : color.match(/\d+/g).slice(0, 3).map(Number));

function stops(points) {
  const list = points.map(([v, color]) => [v, Array.isArray(color) ? color : hex(color)]);
  return (value, rgb) => {
    let n = 0;
    while (n < list.length - 2 && value >= list[n + 1][0]) n++;
    const [v0, a] = list[n], [v1, b] = list[n + 1], t = clamp((value - v0) / (v1 - v0), 0, 1);
    for (let c = 0; c < 3; c++) rgb[c] = (a[c] + t * (b[c] - a[c])) / 255;
    return rgb;
  };
}

const signed = (v, digits) => v.toFixed(digits).replace(/^-(?=.*[1-9])/, '−').replace(/^-/, '');
const dateText = (day) => { const { month, date } = calendarDate(day); return `${month} ${date}`; };

export function mountMovingsea(root) {
  let spinupDays = 3650;
  const runDay = (t) => Math.round(t / 86400 + spinupDays);
  const speedRamp = stops([[0, DARK_NEUTRAL], [0.5, CURRENT], [1, PALE]]), depthRamp = stops([[0, DARK_NEUTRAL], [0.5, VAPOR], [1, PALE]]), rainRamp = stops([[0, DARK_NEUTRAL], [0.5, VAPOR], [1, PALE]]);
  const fields = [
    { name: 'seaSurfaceTemperature', format: signed, label: 'sea surface temperature', unit: '°C', digits: 1, color: (v, rgb) => rampRGB((v + 2) / 34, rgb, DARK_NEUTRAL), legend: [['ramp', 'temperature of the sea’s surface, from −2 to 32\u00a0°C', 'cool', 'warm', DARK]] },
    {
      name: 'current', format: signed, label: 'current', unit: 'm/s', digits: 2, meanDigits: 3,
      compute: (read, out) => { const east = read('currentEast'), north = read('currentNorth'); for (let i = 0; i < out.length; i++) out[i] = Math.hypot(east[i], north[i]); return out; },
      arrow: (read, i, out) => {
        const east = read('currentEast')[i], north = read('currentNorth')[i], speed = Math.hypot(east, north), length = speed ? 0.17 / (speed + 0.25) : 0;
        out[0] = east * length; out[1] = north * length;
      },
      color: speedRamp,
      legend: [['ramp', 'speed of the top layer’s current, from 0 to 1\u00a0m/s', DARK, PALE, CURRENT], ['arrow', 'its direction, the arrows longer for faster currents', 'rgba(255, 255, 255, 0.9)']],
    },
    { name: 'mixedLayerDepth', format: signed, label: 'mixed-layer depth', unit: 'm', digits: 0, color: (v, rgb) => depthRamp(Math.log(Math.max(v, 2) / 2) / Math.log(256), rgb), legend: [['ramp', 'depth of the mixed layer, from 2\u00a0m to 500\u00a0m on a logarithmic scale', DARK, PALE, VAPOR]] },
    { name: 'seaSurfaceHeight', format: signed, label: 'sea level', unit: 'm', digits: 2, mean: false, color: (v, rgb) => rampRGB(0.5 + v / 3, rgb, DARK_NEUTRAL), legend: [['ramp', 'sea surface 1.5\u00a0m below to 1.5\u00a0m above its mean', 'cool', 'warm', DARK]] },
    { name: 'upwelling50', format: signed, label: 'upwelling at 50\u00a0m', unit: 'm a day', digits: 2, mean: false, color: (v, rgb) => rampRGB(0.5 - v / 4, rgb, DARK_NEUTRAL), legend: [['ramp', 'water 50\u00a0m down rising 2\u00a0m a day to sinking 2\u00a0m a day', 'cool', 'warm', DARK]] },
    { name: 'surfaceSalinity', format: signed, label: 'surface salinity', unit: 'psu', digits: 2, meanDigits: 2, color: (v, rgb) => rampRGB((v - 33.6) / 2, rgb, DARK_NEUTRAL), legend: [['ramp', 'salinity of the top layer, from 33.6 to 35.6\u00a0psu', 'cool', 'warm', DARK]] },
    { name: 'rain', format: signed, label: 'rain', unit: 'mm a day', digits: 1, meanDigits: 2, color: (v, rgb) => rainRamp(Math.sqrt(Math.max(0, v) / 30), rgb), legend: [['ramp', 'rain, from none to 30\u00a0mm a day and more, averaged since the frame before', DARK, PALE, VAPOR]] },
  ];
  playbackGlobe(root, {
    url: DATA_URL, N: 16, center: { lat: 15 * Math.PI / 180, lon: 0 },
    sequences: [{ name: 'year', label: 'the year, a frame every five days', framesPerSecond: 6 }, { name: 'spinup', label: 'the first 60 days, a frame a day', framesPerSecond: 6 }],
    fields, globalMean: true, arrows: { coarse: 8, opacity: 0.85 },
    when: {
      slider: (t, sequence) => { const day = t / 86400; return sequence === 'year' ? `day ${Math.round(day)}, ${dateText(day)}` : `day ${runDay(t)} of the run, ${dateText(day)}`; },
      readout: (t, sequence) => {
        const day = t / 86400, { north, south } = season(day);
        return [['day', sequence === 'year' ? `${Math.round(day)}, ${dateText(day)}` : `${runDay(t)} of the run, ${dateText(day)}`], ['season', `${north} in the north, ${south} in the south, ${seasonPhrase(day)}`]];
      },
    },
    note: 'Computed in advance: Part 7’s wet planet, with the climate model’s own winds, radiation, evaporation, convection and clouds on the N = 16 grid with 36 layers, its slab replaced by the model’s own ocean, 58 layers of water over a flat floor 4,000\u00a0m down. There is no land, so the sea has no coasts. It shows one model year from the March equinox after ten years of spin-up, and the first 60 days of the spin-up, a frame a day, as the sea starts moving from rest. Temperature, salinity and current are the top layer’s, 3\u00a0m thick; the mixed-layer depth is the one the ocean measures after it has mixed the water; the rain is an average since the frame before. The run’s starting ice, its bare ice as dark as the sea, caps the polar seas in the first 60 days, but almost all of it melted in the first years. Drag the globe to turn it, or the time slider to step through the days.',
  });
  loadFrames(DATA_URL).then(({ header }) => { spinupDays = header.spinupDays; }).catch(() => {});
}
