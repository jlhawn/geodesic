import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { terrainColor } from "./terrain.module.js";

const SUBDIVISIONS = 320;

export default function runApp() {
  const grid = new Grid(SUBDIVISIONS);
  
  // Initialize Viewer
  initUnifiedViewer(document.body, grid, {
    backgroundColor: 0x222222,
    getColor: (cell) => {
      // if (cell.isPole) {
      //   return {r: 0, g: 1, b: 1};
      // } else if (cell.isPentagon) {
      //   return {r: 1, g: 0, b: 0};
      // } else if (cell.isAlongIcosahedronEdge) {
      //   return {r: 1, g: 1, b: 0};
      // }

      return terrainColor(cell.centerVertex.x, cell.centerVertex.y, cell.centerVertex.z)
    },
  });
};

