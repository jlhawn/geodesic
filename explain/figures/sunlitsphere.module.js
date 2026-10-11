import { playbackGlobe } from '../playbackGlobe.module.js';
import { LAYER_FIELDS } from '../layerCases.module.js';
import { rampRGB, sequentialRGB, termColor, clamp, DARK_NEUTRAL } from '../runtime.module.js';
import { DATA_URL, calendarDate, season } from '../sunlitRun.module.js';
import { seasonPhrase } from '../insolation.module.js';

const DARK = 'rgb(52, 55, 62)', PALE = '#d4efff';

function infraredRamp() {
  const hex = (color) => color.match(/[0-9a-f]{2}/gi).map((h) => parseInt(h, 16));
  const stops = [[200, DARK_NEUTRAL], [290, hex(termColor('b'))], [380, hex(PALE)]];
  return (v, rgb) => {
    const n = v < stops[1][0] ? 0 : 1, [v0, a] = stops[n], [v1, b] = stops[n + 1], t = clamp((v - v0) / (v1 - v0), 0, 1);
    for (let c = 0; c < 3; c++) rgb[c] = (a[c] + t * (b[c] - a[c])) / 255;
    return rgb;
  };
}

const dateText = (day) => { const { month, date } = calendarDate(day); return `${month} ${date}`; };
const hour = (day) => `${String(Math.round(calendarDate(day).hour)).padStart(2, '0')}:00`;
const clock = (day) => `${hour(day)} at Greenwich`;

export function mountSunlitsphere(root) {
  const infrared = termColor('b');
  const fields = [
    { name: 'surfaceTemperature', label: 'ground temperature', unit: 'K', digits: 0, color: (v, rgb) => rampRGB((v - 200) / 100, rgb, DARK_NEUTRAL), legend: [['ramp', 'temperature of the ground, from 200 to 300\u00a0K', 'cool', 'warm', DARK]] },
    { name: 'outgoingLongwave', label: 'infrared leaving the top', unit: 'W/m²', digits: 0, color: infraredRamp(), legend: [['ramp', 'infrared leaving the top of the air, from 200 to 380\u00a0W/m²', DARK, PALE, infrared]] },
    { ...LAYER_FIELDS[1], label: 'temperature 1.5\u00a0km up', color: (v, rgb) => rampRGB((v - 220) / 70, rgb, DARK_NEUTRAL), legend: [['ramp', 'temperature at 850\u00a0hPa, about 1.5\u00a0km up, from 220 to 290\u00a0K', 'cool', 'warm', DARK]] },
    { ...LAYER_FIELDS[2], label: 'wind speed 10\u00a0km up', color: (v, rgb) => sequentialRGB(v / 60, rgb), legend: [['ramp', 'wind speed at 250\u00a0hPa, about 10\u00a0km up, from 0 to 60\u00a0m/s', DARK, 'warm']] },
  ];
  playbackGlobe(root, {
    url: DATA_URL, N: 16, center: { lat: 35 * Math.PI / 180, lon: 0 },
    sequences: [{ name: 'year', label: 'the year, a frame every five days', framesPerSecond: 6 }, { name: 'weather', label: 'a month of weather, four frames a day', framesPerSecond: 8 }],
    fields, globalMean: true,
    when: {
      slider: (t, sequence) => { const day = t / 86400; return sequence === 'year' ? `day ${Math.round(day)}, ${dateText(day)}` : `day ${Math.round(day)}, ${dateText(day)}, ${hour(day)}`; },
      readout: (t, sequence) => {
        const day = t / 86400, { north, south } = season(day);
        return [['day', sequence === 'year' ? `${Math.round(day)}, ${dateText(day)}` : `${Math.round(day)}, ${dateText(day)}, ${clock(day)}`], ['season', `${north} in the north, ${south} in the south, ${seasonPhrase(day)}`]];
      },
    },
    note: 'Computed in advance: Part 5’s dry sphere, 27 layers on the N = 16 grid, heated by the model’s own radiation under the sun of each hour of the year, over a ground that reflects 7 percent of the sunlight. It shows one model year, starting at the March equinox after 150 days with the sun on, and the January within it at four frames a day. Drag the globe to turn it, or the time slider to step through the days.',
  });
}
