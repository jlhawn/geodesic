import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { terrainColor } from "./terrain.module.js";

const SUBDIVISIONS = 100;

export default function runApp() {
  const grid = new Grid(SUBDIVISIONS);
  
  // Initialize Viewer
  initUnifiedViewer(document.body, grid, {
    backgroundColor: 0x222222,
    getColor: (cell) => terrainColor(cell.centerVertex.x, cell.centerVertex.y, cell.centerVertex.z),
  });
};

