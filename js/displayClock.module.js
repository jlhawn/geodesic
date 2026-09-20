/*
 * A display clock that follows the model clock smoothly. The model
 * advances in jumps of many minutes a few times a second; the display
 * clock advances every animation frame at the measured simulation rate
 * and aims to stay `lead` frame intervals behind the latest model time.
 * A low-gain proportional controller nudges its speed: a little faster
 * when it has fallen further behind, a little slower as it catches up,
 * so a frame's arrival changes the speed by a few percent rather than
 * doubling it. It never runs more than `lead` intervals ahead of the
 * model, snaps after a jump of many frames, and while the model is
 * paused it eases forward onto the model time without ever running
 * backwards.
 */
export function createDisplayClock({ lead = 2, gain = 0.1, snapFrames = 20 } = {}) {
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
          const factor = Math.min(2, Math.max(0, 1 + gain * (error - lag) / lag)) * Math.min(1, Math.max(0, (error + lag) / lag));
          time += dt * rate * factor;
        }
      } else if (error > 0) {
        time += error * Math.min(1, dt / 500);
      }
      return time;
    },
  };
}
