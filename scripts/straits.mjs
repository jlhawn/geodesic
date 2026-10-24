/*
 * Passages between basins for the strait checks: a point inside each
 * basin (latitude, longitude) and a box [south, north, west, east] that
 * holds the passage and no other route between them.
 */
export const STRAITS = {
  gibraltar: { a: [36.0, -9.0], b: [37.0, -1.0], box: [34.5, 38.5, -10, 0] },
  babElMandeb: { a: [16.5, 40.5], b: [12.0, 45.0], box: [11, 17, 40, 46] },
  hormuz: { a: [26.6, 53.5], b: [25.0, 57.8], box: [23.8, 27.8, 52.5, 59] },
  bering: { a: [62.5, -172], b: [67.5, -168.5], box: [62, 69, -177, -163] },
  lombok: { a: [-8.2, 115.85], b: [-9.0, 115.8], box: [-9.6, -7.8, 115.5, 116.1] },
  panama: { a: [8.0, -79.5], b: [9.7, -79.5], box: [6.5, 11, -84, -76] },
  danish: { a: [57.3, 11.5], b: [55.3, 15.5], box: [54.2, 58.0, 9.5, 16.5] },
  bosporus: { a: [42.2, 29.5], b: [40.8, 28.3], box: [40.4, 42.6, 27.5, 30.5] },
  hudson: { a: [60.0, -84.0], b: [61.0, -62.5], box: [57.5, 65.5, -86, -61] },
  tsushima: { a: [37.0, 131.0], b: [32.5, 127.5], box: [31.5, 37.5, 126.5, 131.3] },
  tatar: { a: [47.0, 139.8], b: [54.2, 141.2], box: [46.5, 55, 137, 142.3] },
  florida: { a: [25.0, -85.0], b: [27.5, -79.8], box: [22.6, 28.5, -86, -79.0] },
  denmarkStrait: { a: [64.5, -31.0], b: [68.0, -22.0], box: [63.5, 69.0, -36, -20.5] },
  faroeBank: { a: [61.8, -3.0], b: [61.0, -12.0], box: [59.5, 62.4, -13, -2.5] },
};
