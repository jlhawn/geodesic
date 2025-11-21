
const EPS = 1e-6;

export function lerp(aVal, aMin, aMax, bMin, bMax) {
    const aRange = aMax - aMin;
    const valRatio = (aVal - aMin)/aRange;

    const bRange = bMax - bMin;
    const bVal = valRatio*bRange + bMin;
    return bVal;
};

export function snapNearZero(x) {
    return Math.abs(x) < EPS ? 0 : x;
}

export function clamp(val, min, max) {
    return Math.min(max, Math.max(min, val));
};
