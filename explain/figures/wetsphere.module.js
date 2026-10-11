import { playbackGlobe } from '../playbackGlobe.module.js';
import { rampRGB, termColor, clamp, DARK_NEUTRAL, VAPOR, CLOUD } from '../runtime.module.js';
import { DATA_URL, calendarDate, season } from '../wetRun.module.js';
import { seasonPhrase } from '../insolation.module.js';

const DARK = 'rgb(52, 55, 62)', PALE = '#d4efff';
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

const dateText = (day) => { const { month, date } = calendarDate(day); return `${month} ${date}`; };
const hour = (day) => `${String(Math.round(calendarDate(day).hour)).padStart(2, '0')}:00`;
const clock = (day) => `${hour(day)} at Greenwich`;

export function mountWetsphere(root) {
  const infrared = termColor('b');
  const rainRamp = stops([[0, DARK_NEUTRAL], [0.5, VAPOR], [1, PALE]]), cloudRamp = stops([[0, DARK_NEUTRAL], [1, CLOUD]]);
  const fields = [
    { name: 'rain', label: 'rain', unit: 'mm a day', digits: 1, meanDigits: 2, color: (v, rgb) => rainRamp(Math.sqrt(Math.max(0, v) / 30), rgb), legend: [['ramp', 'rain, from none to 30\u00a0mm a day and more, brighter for heavier', DARK, PALE, VAPOR]] },
    { name: 'precipitableWater', label: 'water in the air', unit: 'kg/m²', digits: 0, color: stops([[0, DARK_NEUTRAL], [35, VAPOR], [60, PALE]]), legend: [['ramp', 'vapor over each square meter, from 0 to 60\u00a0kg/m²', DARK, PALE, VAPOR]] },
    { name: 'cloudWater', label: 'cloud water', unit: 'g/m²', digits: 0, color: (v, rgb) => cloudRamp(Math.sqrt(Math.max(0, v) / 400), rgb), legend: [['ramp', 'cloud water over each square meter, from 0 to 400\u00a0g/m² and more', DARK, CLOUD]] },
    { name: 'surfaceTemperature', label: 'surface temperature', unit: 'K', digits: 0, color: (v, rgb) => rampRGB((v - 250) / 55, rgb, DARK_NEUTRAL), legend: [['ramp', 'temperature of the sea or its ice, from 250 to 305\u00a0K', 'cool', 'warm', DARK]] },
    { name: 'outgoingLongwave', label: 'infrared leaving the top', unit: 'W/m²', digits: 0, color: stops([[100, DARK_NEUTRAL], [230, infrared], [320, PALE]]), legend: [['ramp', 'infrared leaving the top of the air, from 100 to 320\u00a0W/m²', DARK, PALE, infrared]] },
  ];
  playbackGlobe(root, {
    url: DATA_URL, N: 16, center: { lat: 0, lon: 0 },
    sequences: [{ name: 'year', label: 'the year, a frame every five days', framesPerSecond: 6 }, { name: 'weather', label: 'a month of weather, four frames a day', framesPerSecond: 8 }],
    fields, globalMean: true,
    when: {
      slider: (t, sequence) => { const day = t / 86400; return sequence === 'year' ? `day ${Math.round(day)}, ${dateText(day)}` : `day ${Math.round(day)}, ${dateText(day)}, ${hour(day)}`; },
      readout: (t, sequence) => {
        const day = t / 86400, { north, south } = season(day);
        return [['day', sequence === 'year' ? `${Math.round(day)}, ${dateText(day)}` : `${Math.round(day)}, ${dateText(day)}, ${clock(day)}`], ['season', `${north} in the north, ${south} in the south, ${seasonPhrase(day)}`]];
      },
    },
    note: 'Computed in advance: the climate model’s own winds, radiation, evaporation, convection and clouds on the N = 16 grid with 36 layers, over a planet of open sea with no land and no currents, a slab of water 10\u00a0m deep that freezes into sea ice, its bare ice made as dark as the sea. It shows one model year, starting at the March equinox after a year of spin-up, and the January within it at four frames a day. The rain of each frame is its average since the frame before. Drag the globe to turn it, or the time slider to step through the days.',
  });
}
