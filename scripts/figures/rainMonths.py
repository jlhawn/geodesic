# Three-row rain figure of a quarter's saves from rainMonths.mjs: each save's 30-day rain total in mm
# beside the model minus GPCP over the same window.  python3 scripts/figures/rainMonths.py <json> <png> [title]
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from plotting import load, default_title, array, Raster, style_axes, world_ticks, BACKGROUND, TEXT, MUTED

d = load(sys.argv[1])
title = sys.argv[3] if len(sys.argv) > 3 else f"{d['tag']} N={d['N']}, day {d['day']}: rain by save against GPCP"
lon, lat, land = array(d['lon']), array(d['lat']), array(d['land']) > 0.5
raster = Raster(lon, lat, 0.25 if d['N'] > 64 else 0.5)
landG = raster(land.astype(float))

RAIN_BOUNDS = [0, 10, 25, 50, 100, 200, 400, 800]
rain_cmap = ListedColormap(LinearSegmentedColormap.from_list('rain', ['#101a36', '#16407a', '#1f78b4', '#4fb3e8', '#a6e3ff', '#f4fbff'])(np.linspace(0, 1, len(RAIN_BOUNDS))))
rain_norm = BoundaryNorm(RAIN_BOUNDS, rain_cmap.N, extend='max')

DIFF_BOUNDS = [-200, -100, -50, -25, -10, 10, 25, 50, 100, 200]
brbg = plt.get_cmap('BrBG')
diff_colors = [brbg(x) for x in np.linspace(0.08, 0.92, len(DIFF_BOUNDS) + 1)]
diff_colors[len(DIFF_BOUNDS) // 2] = '#2a3042'
diff_cmap = ListedColormap(diff_colors)
diff_norm = BoundaryNorm(DIFF_BOUNDS, diff_cmap.N, extend='both')

fig, axes = plt.subplots(len(d['saves']), 2, figsize=(16, 4.1 * len(d['saves']) + 1.6), facecolor=BACKGROUND, constrained_layout=True)
mm = lambda x: f'{x:.0f}'
sg = lambda x: f'{x:+.0f}'
for row, s in zip(axes, d['saves']):
    model, obs = array(s['model']), array(s['obs'])
    st = s['stats']
    for ax, field, cmap, norm in ((row[0], model, rain_cmap, rain_norm), (row[1], model - obs, diff_cmap, diff_norm)):
        style_axes(ax); world_ticks(ax)
        im = ax.imshow(np.ma.masked_invalid(raster(field)), extent=raster.extent, origin='lower', cmap=cmap, norm=norm, interpolation='nearest')
        ax.contour(raster.gx, raster.gy, landG, levels=[0.5], colors='#c8ccd6', linewidths=0.35)
        if ax is row[0]: left = im
        else: right = im
    g, l, o = st['global'], st['land'], st['sea']
    row[0].set_title(f"Day {s['day']}, {s['window']} ({s['weights']})\nmodel rain: global {mm(g['model'])} mm ({g['model'] / s['days']:.2f} mm/d), land {mm(l['model'])}, sea {mm(o['model'])}", color=TEXT, fontsize=10, loc='left')
    row[1].set_title(f"Day {s['day']}, model minus GPCP: global {sg(g['diff'])} mm, land {sg(l['diff'])}, sea {sg(o['diff'])}; mean |difference| {mm(g['absDiff'])}\nGPCP: global {mm(g['obs'])} mm, land {mm(l['obs'])}, sea {mm(o['obs'])}", color=TEXT, fontsize=10, loc='left')
for im, ax, label, ticks in ((left, axes[:, 0], 'Rain over the 30-day window, mm', RAIN_BOUNDS), (right, axes[:, 1], 'Model minus GPCP over the same window, mm (brown too dry, teal too wet)', DIFF_BOUNDS)):
    cb = fig.colorbar(im, ax=ax, location='bottom', shrink=0.8, aspect=45, pad=0.02, ticks=ticks, spacing='uniform')
    cb.ax.xaxis.set_tick_params(color=MUTED, labelsize=8)
    plt.setp(cb.ax.get_xticklabels(), color=MUTED)
    cb.outline.set_edgecolor('#444')
    cb.set_label(label, color=TEXT, fontsize=9)
fig.suptitle(title, color=TEXT, fontsize=13)
fig.supxlabel('GPCP v2.3 is a 1991–2020 monthly climatology weighted by the days of each calendar month in the window; each model window is a single model year.', color=MUTED, fontsize=8.5)
fig.savefig(sys.argv[2], dpi=100 if d['N'] > 64 else 80, facecolor=fig.get_facecolor())
print('wrote', sys.argv[2])
