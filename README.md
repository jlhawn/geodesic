# WebGCM

A global climate model that runs in the browser. The repository keeps the name of the geodesic grid the model is built on.

Live at https://gcm.echorelay.net/

- **[Climate model](https://gcm.echorelay.net/climate.html)** — a hydrostatic atmosphere on a geodesic C-grid with moist physics, a six-layer ocean, sea ice, snow and real continents with terrain, integrated on the GPU through WebGPU. Design notes in [docs/c-grid-dynamical-core.md](docs/c-grid-dynamical-core.md), plans in [docs/roadmap.md](docs/roadmap.md).
- **[Geodesic grid](https://gcm.echorelay.net/grid.html)** — the original geodesic polyhedron viewer (Three.js).
- **[Synoptic charts](https://gcm.echorelay.net/charts.html)** — pressure-level charts drawn from saved model states.

To run locally, serve the repository with `python3 httpd.py` and open http://localhost:8000/. The model needs the page cross-origin isolated (`Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`), which that server and the `_headers` file for Cloudflare or Netlify Pages both provide. `npm test` runs the test suite.

![The simulated Earth from space](screenshots/satellite.jpg)
