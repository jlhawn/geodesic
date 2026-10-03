# Four global maps of a state from stateMaps.mjs: mixed-layer temperature, sea-ice thickness, vegetation cover and
# surface temperature.  python3 scripts/figures/stateMaps.py <json> <png> [title]
import sys
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from plotting import load, default_title, array, Raster, style_axes, world_ticks, colourbar, BACKGROUND, LAND, TEXT

d = load(sys.argv[1])
title = sys.argv[3] if len(sys.argv) > 3 else default_title(d)
lon, lat = array(d['lon']), array(d['lat'])
land = array(d['land']) > 0.5
ice, conc, snow = array(d['ice']), array(d['concentration']), array(d['snow'])
sst, ts, veg, trees = array(d['sst']), array(d['ts']), array(d['vegetation']), array(d['trees'])
raster = Raster(lon, lat, 0.25 if d['N'] > 64 else 0.5)
landG = raster(land.astype(float)) > 0.5
iced = ~land & (ice > 0)
fig, axes = plt.subplots(2, 2, figsize=(16, 8.6), facecolor=BACKGROUND)


def panel(ax, field, cmap, vmin, vmax, label, under, extend='neither', background=None):
    style_axes(ax); world_ticks(ax)
    if background is not None:
        ax.imshow(background, extent=raster.extent, origin='lower', cmap=ListedColormap(['#05070d', LAND]), vmin=0, vmax=1, interpolation='nearest')
    im = ax.imshow(np.ma.masked_invalid(field), extent=raster.extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    colourbar(fig, im, ax, under, fraction=0.025, pad=0.01, extend=extend)
    ax.set_title(label, color=TEXT, fontsize=10.5, loc='left')


seaBackground = landG.astype(float)
panel(axes[0, 0], np.where(landG, np.nan, raster(sst)), 'RdYlBu_r', -2, 30, 'Mixed-layer temperature °C', '°C', extend='both', background=seaBackground)
iceG = raster(np.where(iced, ice, np.nan))
meanSnow = float(np.mean(snow[iced])) if iced.any() else 0.0
panel(axes[0, 1], np.where(landG, np.nan, iceG), 'Blues', 0, 3, f"Sea-ice thickness m (iced cells {int(iced.sum())}, mean snow {meanSnow:.0f} kg/m² over them)", 'm', extend='max', background=seaBackground)
landMean = lambda v: float(np.mean(v[land])) if land.any() else 0.0
panel(axes[1, 0], np.where(landG, raster(veg), np.nan), 'YlGn', 0, 1, f"Vegetation cover v, 0 bare – 1 closed (trees or grass) (land means: v {landMean(veg):.2f}, trees {landMean(trees):.2f})", 'cover', background=seaBackground)
panel(axes[1, 1], raster(ts), 'coolwarm', -30, 40, 'Surface temperature °C', '°C', extend='both')
fig.suptitle(title, color=TEXT, fontsize=13)
plt.tight_layout(rect=(0, 0, 1, 0.965))
fig.savefig(sys.argv[2], dpi=100 if d['N'] > 64 else 80, facecolor=fig.get_facecolor())
print('wrote', sys.argv[2])
