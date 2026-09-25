/*
 * How a pair of touch points moved between two samples, in screen
 * pixels: the shift of their midpoint (dx, dy), the ratio of their
 * spread (scale), the turn of the line through them in radians,
 * clockwise on screen and within ±π (turn), and the new midpoint (x, y).
 */
export function twoFingerMotion(a0, b0, a1, b1) {
  const x = (a1.x + b1.x) / 2, y = (a1.y + b1.y) / 2;
  const spread0 = Math.hypot(b0.x - a0.x, b0.y - a0.y), spread1 = Math.hypot(b1.x - a1.x, b1.y - a1.y);
  let turn = Math.atan2(b1.y - a1.y, b1.x - a1.x) - Math.atan2(b0.y - a0.y, b0.x - a0.x);
  if (turn > Math.PI) turn -= 2 * Math.PI;
  else if (turn < -Math.PI) turn += 2 * Math.PI;
  return { dx: x - (a0.x + b0.x) / 2, dy: y - (a0.y + b0.y) / 2, scale: spread0 > 0 && spread1 > 0 ? spread1 / spread0 : 1, turn, x, y };
}
