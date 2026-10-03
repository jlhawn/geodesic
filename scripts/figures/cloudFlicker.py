# The blink figures of scripts/cloudFlicker.mjs: a map of how often each cell's cloud overlay blinks, and a strip of
# eight consecutive steps and eight of the page's frames over the box where blinks are densest, painted as the page's
# Cloud cover overlay paints them (white at opacity 1 - exp(-g / 40 g/m^2) over sea blue and land brown).
#   python3 scripts/figures/cloudFlicker.py <plot.json> <map.png> <strip.png>
import sys
import numpy as np
import matplotlib.pyplot as plt
from plotting import load, default_title, array, Raster, style_axes, world_ticks, colourbar, BACKGROUND, TEXT, MUTED

d = load(sys.argv[1])
title = default_title(d)
lon, lat = array(d['lon']), array(d['lat'])
land = array(d['land']) > 0.5
frequency = array(d['blinks']) / max(1, d['transitions'])
hours = d['dt'] / 3600

raster = Raster(lon, lat, 0.25 if d['N'] > 64 else 0.5)
fig, ax = plt.subplots(figsize=(15, 7.6), facecolor=BACKGROUND)
style_axes(ax); world_ticks(ax)
F = raster(np.where(frequency > 0, frequency, np.nan))
ax.imshow(raster(land.astype(float)), extent=raster.extent, origin='lower', cmap=plt.matplotlib.colors.ListedColormap(['#10182b', '#3a3a3a']), vmin=0, vmax=1, interpolation='nearest')
im = ax.imshow(np.ma.masked_invalid(F), extent=raster.extent, origin='lower', cmap='inferno', vmin=0, vmax=max(0.05, float(np.nanpercentile(F, 99.5)) if np.isfinite(F).any() else 0.05), interpolation='nearest')
colourbar(fig, im, ax, 'blink onsets per step', fraction=0.025, pad=0.01, extend='max')
s, n, w, e = d['strip']['box']
for x0, x1 in ([(w, e)] if w <= e else [(w, 180), (-180, e)]):
    ax.plot([x0, x1, x1, x0, x0], [s, s, n, n, s], color='#7fd1ff', lw=1)
for x, y in d['worst']:
    ax.plot(x, y, marker='o', mfc='none', mec='#7fd1ff', ms=7, mew=1)
blinking = float(np.mean(frequency > 0))
ax.set_title(f"Cloud-overlay blinks: opacity moves > {d['jump']} in one step and back within three ({d['transitions']} counted steps of {d['dt']:.2f} s; "
             f"{100 * blinking:.1f}% of cells blink at least once; mean {100 * float(np.mean(frequency)):.2f}% of cells per step)",
             color=TEXT, fontsize=10, loc='left')
fig.suptitle(f"{title}: blink frequency (box: strip region, rings: the worst cells)", color=TEXT, fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.965))
fig.savefig(sys.argv[2], dpi=100 if d['N'] > 64 else 80, facecolor=fig.get_facecolor())

st = d['strip']
centre = (w + (e if w <= e else e + 360)) / 2
slon = (array(st['lon']) - centre + 180) % 360 - 180
half = ((e if w <= e else e + 360) - w) / 2
box = Raster(slon, array(st['lat']), 0.1 if d['N'] > 64 else 0.2, window=(-half, half, s, n))
landBox = box(array(st['land'])) > 0.5
OCEAN = np.array([0.05, 0.22, 0.45]); LANDC = np.array([0.30, 0.33, 0.17]); WHITE = np.ones(3)


def painted(alpha):
    A = np.nan_to_num(box(array(alpha)))[..., None]
    base = np.where(landBox[..., None], LANDC, OCEAN)
    return base * (1 - A) + WHITE * A


rows = [('step', st['frames'], [f"+{h:.2f} h" for h in st['hours']]),
        ('page frame', st['pageFrames'], [f"+{(d['stride'] * (k + 1)) * hours:.2f} h" for k in range(len(st['pageFrames']))])]
fig, axes = plt.subplots(2, 8, figsize=(22, 5.2), facecolor=BACKGROUND)
for r, (label, frames, stamps) in enumerate(rows):
    for c in range(8):
        a = axes[r, c]
        style_axes(a)
        a.set_xticks([]); a.set_yticks([])
        if c >= len(frames):
            a.set_visible(False); continue
        a.imshow(painted(frames[c]), extent=box.extent, origin='lower', interpolation='nearest')
        a.set_title(f"{label} {stamps[c]}", color=MUTED, fontsize=8.5)
lonLabel = lambda x: f"{abs(x):.0f}{'E' if x >= 0 else 'W'}"
fig.suptitle(f"{title}: the Cloud cover overlay over {abs(s):.0f}{'S' if s < 0 else 'N'}–{abs(n):.0f}{'S' if n < 0 else 'N'}, {lonLabel(w)}–{lonLabel(e)}; "
             f"top: {len(st['frames'])} consecutive steps ({d['dt']:.2f} s each); bottom: one step in {d['stride']}, as the page shows its frames",
             color=TEXT, fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(sys.argv[3], dpi=100, facecolor=fig.get_facecolor())
