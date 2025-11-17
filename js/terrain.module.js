import SimplexNoise from "./simplex-noise.module.js";

const simplex = new SimplexNoise(Math.random());

/* ===== TUNING CONSTANTS ===== */

// Continents (land vs ocean)
const CONT_SCALE        = 1.4;   // frequency of continent field (1.2–1.6 is a good range)
const CONT_OCTAVES      = 2;     // fewer octaves → simpler, blob-like continents
const CONT_LACUNARITY   = 1.9;
const CONT_PERSISTENCE  = 0.6;

// Coastline “wiggle” (adds fine detail to coasts, but should not reconnect continents)
const COAST_NOISE_SCALE      = 4.0;
const COAST_NOISE_OCTAVES    = 1;
const COAST_NOISE_LACUNARITY = 2.0;
const COAST_NOISE_PERSIST    = 0.5;
const COAST_NOISE_AMPLITUDE  = 0.03; // small; raise slightly for more jagged coasts

// Sea level: higher → more ocean, smaller continents
const SEA_LEVEL = 0.55; // 0.6–0.7: more ocean, 0.5–0.6: more land

// Ocean depth scaling
const OCEAN_DEPTH_SCALE = 1.0;

// Land detail – low and high frequency
const DETAIL_LOW_SCALE       = 2.0;
const DETAIL_LOW_OCTAVES     = 3;
const DETAIL_LOW_LACUNARITY  = 2.0;
const DETAIL_LOW_PERSISTENCE = 0.5;

const DETAIL_HIGH_SCALE       = 6.0;
const DETAIL_HIGH_OCTAVES     = 2;
const DETAIL_HIGH_LACUNARITY  = 2.5;
const DETAIL_HIGH_PERSISTENCE = 0.5;

// Mountain ridges
const RIDGED_SCALE       = 5.0;
const RIDGED_OCTAVES     = 3;
const RIDGED_LACUNARITY  = 2.0;
const RIDGED_PERSISTENCE = 0.5;

// Land shaping
const LAND_SMOOTHSTEP = true;   // smooth landMask edges
const PLAINS_FACTOR   = 0.9;    // base land contribution
const DETAIL_LOW_WEIGHT  = 0.25;
const DETAIL_HIGH_WEIGHT = 0.18;
const MOUNTAIN_WEIGHT    = 0.35;
const LAND_EXPONENT      = 1.2; // >1 makes high peaks rarer


/* ===== NOISE HELPERS ===== */

function fbmSimplex(x, y, z, octaves, lacunarity, persistence) {
  let frequency = 1.0;
  let amplitude = 1.0;
  let value = 0.0;
  let maxAmp = 0.0;

  for (let i = 0; i < octaves; i++) {
    value += simplex.noise3D(x * frequency, y * frequency, z * frequency) * amplitude;
    maxAmp += amplitude;

    frequency *= lacunarity;
    amplitude *= persistence;
  }

  // Normalize to roughly [-1, 1]
  return value / maxAmp;
}

function ridgedSimplex(x, y, z, octaves, lacunarity, persistence) {
  let frequency = 1.0;
  let amplitude = 1.0;
  let value = 0.0;
  let maxAmp = 0.0;

  for (let i = 0; i < octaves; i++) {
    let n = simplex.noise3D(x * frequency, y * frequency, z * frequency);
    n = 1.0 - Math.abs(n);   // ridged: peaks along noise creases
    value += n * amplitude;
    maxAmp += amplitude;

    frequency *= lacunarity;
    amplitude *= persistence;
  }

  // Roughly [0, 1]
  return value / maxAmp;
}


/* ===== MAIN TERRAIN FUNCTION ===== */

function earthLikeHeight(x, y, z) {
  // 1. Continents: low-frequency mask → decides land vs ocean.
  let cont = fbmSimplex(
    x * CONT_SCALE,
    y * CONT_SCALE,
    z * CONT_SCALE,
    CONT_OCTAVES,
    CONT_LACUNARITY,
    CONT_PERSISTENCE
  ); // ~[-1, 1]

  // Map to [0, 1]
  let cont01 = (cont + 1) * 0.5;

  // 1b. Coastline “wiggle” – small high-frequency perturbation
  if (COAST_NOISE_AMPLITUDE !== 0) {
    const coastNoise = fbmSimplex(
      x * COAST_NOISE_SCALE,
      y * COAST_NOISE_SCALE,
      z * COAST_NOISE_SCALE,
      COAST_NOISE_OCTAVES,
      COAST_NOISE_LACUNARITY,
      COAST_NOISE_PERSIST
    ); // ~[-1, 1]

    cont01 += COAST_NOISE_AMPLITUDE * coastNoise;
    cont01 = Math.min(1, Math.max(0, cont01)); // clamp
  }

  // 2. Land/ocean split
  const base = cont01 - SEA_LEVEL;     // negative = ocean, positive = land candidate

  const oceanDepth = Math.min(0, base) * OCEAN_DEPTH_SCALE;

  // Land mask [0,1]: 0 at sea level, 1 at highest continents
  let landMask = Math.max(0, base) / (1 - SEA_LEVEL);

  if (LAND_SMOOTHSTEP) {
    // Smooth edges: landMask = smoothstep(0,1,landMask)
    landMask = landMask * landMask * (3 - 2 * landMask);
  }

  // 3. Interior detail on land
  const detailLow = fbmSimplex(
    x * DETAIL_LOW_SCALE,
    y * DETAIL_LOW_SCALE,
    z * DETAIL_LOW_SCALE,
    DETAIL_LOW_OCTAVES,
    DETAIL_LOW_LACUNARITY,
    DETAIL_LOW_PERSISTENCE
  ); // ~[-1, 1]

  const detailHigh = fbmSimplex(
    x * DETAIL_HIGH_SCALE,
    y * DETAIL_HIGH_SCALE,
    z * DETAIL_HIGH_SCALE,
    DETAIL_HIGH_OCTAVES,
    DETAIL_HIGH_LACUNARITY,
    DETAIL_HIGH_PERSISTENCE
  ); // ~[-1, 1]

  // 4. Mountain ridges
  const ridged = ridgedSimplex(
    x * RIDGED_SCALE,
    y * RIDGED_SCALE,
    z * RIDGED_SCALE,
    RIDGED_OCTAVES,
    RIDGED_LACUNARITY,
    RIDGED_PERSISTENCE
  ); // ~[0, 1]

  // 5. Shape land elevation
  let elevation = 0.0;

  // Base land from continent strength
  elevation += landMask * PLAINS_FACTOR;

  // Add low- and high-frequency detail (only where landMask > 0)
  elevation += landMask * DETAIL_LOW_WEIGHT * detailLow;
  elevation += landMask * DETAIL_HIGH_WEIGHT * detailHigh;

  // Mountains from ridged (centered around zero)
  elevation += landMask * MOUNTAIN_WEIGHT * (ridged - 0.5);

  // Clamp land to >= 0, then apply mild non-linearity
  elevation = Math.max(0, elevation);
  elevation = Math.pow(elevation, LAND_EXPONENT);

  // 6. Combine ocean (negative) and land (positive)
  const final = oceanDepth + elevation; // <0 ocean, >0 land, ≈0 coast

  return final;
}

export { earthLikeHeight };
