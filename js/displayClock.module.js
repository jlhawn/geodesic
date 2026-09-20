/*
 * A display clock that follows the model clock smoothly. The model
 * advances in jumps of many minutes a few times a second; the display
 * clock advances every animation frame at the rate the caller measures
 * over the last few model frames, and aims to stay `lead` frame
 * intervals behind the latest model time. Only a small correction,
 * within ±maxDelta of that rate, pulls it back toward the setpoint when
 * it drifts, so the arrival of a frame barely changes its speed. It
 * never runs more than `lead` intervals ahead of a stalled model, snaps
 * after a jump of many frames, and while the model is paused it eases
 * forward onto the model time without ever running backwards.
 */
export function createDisplayClock({ lead = 2, gain = 0.05, maxDelta = 0.1, snapFrames = 20 } = {}) {
  let time = null, wall = null;
  return {
    get time() { return time; },
    advance(now, target, { rate = 0, interval = 500, running = true } = {}) {
      const dt = wall === null ? 0 : now - wall;
      wall = now;
      if (time === null) time = target;
      const error = target - time;
      const lag = lead * interval * rate;
      if (running && rate > 0 && lag > 0) {
        if (Math.abs(error) > snapFrames * interval * rate) time = target;
        else {
          const correction = Math.min(maxDelta, Math.max(-maxDelta, gain * (error - lag) / lag));
          const stall = Math.min(1, Math.max(0, (error + lag) / lag));
          time += dt * rate * (1 + correction) * stall;
        }
      } else if (error > 0) {
        time += error * Math.min(1, dt / 500);
      }
      return time;
    },
  };
}
