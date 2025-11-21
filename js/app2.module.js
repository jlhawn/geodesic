import * as THREE from "./three.module.js";
import { Grid } from "./grid.module.js";
import { initEqualEarthViewer } from "./equalEarthViewer.module.js";
import { terrainColor } from "./terrain.module.js";

// The equation for the number of cells is 10*N^2+2
//  10 ->     1,002
//  20 ->     4,002
//  50 ->    25,002
// 100 ->   100,002
// 200 ->   400,002
// 320 -> 1,024,002
// 400 -> 1,600,002
// 500 -> 2,500,002
const SUBDIVISIONS = 100;
const BACKGROUND_COLOR = 0x222222;

if (SUBDIVISIONS < 2) {
  throw "SUBDIVISIONS must be greater than 2.";
}

export default function runApp() {
    const grid = new Grid(SUBDIVISIONS);

    initEqualEarthViewer(document.body, grid, {
        backgroundColor: BACKGROUND_COLOR,
        getColor: (cell) => {
            const cv = cell.centerVertex;
            return terrainColor(cv.x, cv.y, cv.z);
        },
    });
};
