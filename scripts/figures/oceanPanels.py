# Four panels of the ocean from oceanPanels.mjs: the sea surface temperature,
# its difference from WOA23's monthly climatology (the state's window-mean
# SST against the climatology over the window's days when the dump carries
# one, else the instant against the state's date), the current over the
# top 15 m and the mixed-layer depth.
#   python3 scripts/figures/oceanPanels.py <json> <png> [title]
import sys, numpy as np
import matplotlib; matplotlib.use('Agg')
from plotting import load, default_title, array, Raster, style_axes, world_ticks, colourbar, BACKGROUND, LAND, TEXT, MUTED
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

d = load(sys.argv[1])
title = sys.argv[3] if len(sys.argv) > 3 else default_title(d)
lon, lat = array(d['lon']), array(d['lat'])
land = array(d['land']) > 0.5
sea = ~land
sst, woa, ue, vn, mld, ice = (array(d[k]) for k in ('sst', 'woa', 'ue', 'vn', 'mld', 'ice'))
window = d.get('sstWindow') or {'kind': 'instant'}
windowed = window['kind'] == 'mean'
if windowed:
    first, last = window['days']
    judged, judgedLabel, judgedName = array(d['sstMean']), f"{window['length']}-day mean SST − WOA23, days {first}–{last}", f"{window['length']}-day mean SST of days {first}–{last}"
else:
    judged, judgedLabel, judgedName = sst, 'SST − WOA23 on this date, instant', 'instant SST'
bias = judged - woa
speed = np.hypot(ue, vn)
coslat = np.cos(np.radians(lat))
raster = Raster(lon, lat, 0.25 if d['N'] > 64 else 0.5)
landG = raster(land.astype(float)) > 0.5
seaBackground = landG.astype(float)


def mean(values, mask):
    m = mask & np.isfinite(values)
    return float(np.sum(values[m] * coslat[m]) / np.sum(coslat[m])) if m.any() else float('nan')


def box(lat0, lat1, lon0, lon1):
    w = np.where(lon > 180, lon - 360, lon)
    return sea & (lat >= lat0) & (lat <= lat1) & (w >= lon0) & (w <= lon1)


globalBias = mean(bias, sea & (np.abs(lat) <= 60))
tongue = mean(bias, box(-5, 5, -130, -90))
tropics = mean(bias, sea & (np.abs(lat) <= 20))
fig, axes = plt.subplots(2, 2, figsize=(16, 8.6), facecolor=BACKGROUND)


def panel(ax, field, cmap, vmin, vmax, label, under, extend='neither'):
    style_axes(ax); world_ticks(ax)
    ax.imshow(seaBackground, extent=raster.extent, origin='lower', cmap=ListedColormap(['#05070d', LAND]), vmin=0, vmax=1, interpolation='nearest')
    im = ax.imshow(np.ma.masked_invalid(np.where(landG, np.nan, field)), extent=raster.extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    colourbar(fig, im, ax, under, fraction=0.025, pad=0.01, extend=extend)
    ax.set_title(label, color=TEXT, fontsize=10.5, loc='left')
    return ax


iceG = raster(np.where(sea & (ice > 0), 1.0, np.nan))
panel(axes[0, 0], raster(sst), 'RdYlBu_r', -2, 30, f"Sea surface temperature °C{', instant' if windowed else ''} (ice-covered cells stippled)", '°C', extend='both')
axes[0, 0].contourf(raster.LON, raster.LAT, np.ma.masked_invalid(iceG), levels=[0.5, 1.5], colors='none', hatches=['....'])
panel(axes[0, 1], raster(bias), 'RdBu_r', -5, 5, f"{judgedLabel}, K\n60S–60N {globalBias:+.2f}, 20S–20N {tropics:+.2f}, cold tongue 5S–5N 130–90W {tongue:+.2f}", 'K', extend='both')
ax = panel(axes[1, 0], raster(speed), 'magma', 0, 1.2, f"Current over the top {d.get('surfaceDepth', 15)} m, m/s (max {np.nanmax(speed):.2f})", 'm/s', extend='max')
step = 12 if d['N'] > 64 else 6
U, V = raster(ue), raster(vn)
sl = (slice(step // 2, None, step), slice(step // 2, None, step))
ax.quiver(raster.LON[sl], raster.LAT[sl], U[sl], V[sl], color='#d8dde8', scale=0.08, scale_units='xy', width=0.0012, headwidth=4, headlength=5, minshaft=1.5, alpha=0.85)
panel(axes[1, 1], raster(mld), 'viridis_r', 0, 300, f"Mixed-layer depth m, density 0.03 kg/m³ over the surface layer's (60S–60N mean {mean(mld, sea & (np.abs(lat) <= 60)):.0f}, 20S–20N {mean(mld, sea & (np.abs(lat) <= 20)):.0f})", 'm', extend='max')
fig.suptitle(f"{title}, {d.get('vertical', 'slab')} ocean", color=TEXT, fontsize=13)
plt.tight_layout(rect=(0, 0, 1, 0.965))
fig.savefig(sys.argv[2], dpi=100 if d['N'] > 64 else 80, facecolor=fig.get_facecolor())
print(f"ocean: {judgedName} − WOA 60S–60N {globalBias:+.3f} K, 20S–20N {tropics:+.3f} K, cold tongue {tongue:+.3f} K; max surface current {np.nanmax(speed):.2f} m/s; mixed layer 60S–60N {mean(mld, sea & (np.abs(lat) <= 60)):.1f} m")
