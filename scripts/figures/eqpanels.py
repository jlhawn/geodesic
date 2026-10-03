# Two stacked panels of the equatorial Pacific from eqpanels.mjs, the atmosphere over the ocean.
#   python3 scripts/figures/eqpanels.py <json> <png> [title]
import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
from plotting import default_title
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm
d = json.load(open(sys.argv[1])); cols = d['columns']; A = np.array([[np.nan if v is None else v for v in r] for r in d['rows']], dtype=float)
c = {k: A[:, i] for i, k in enumerate(cols)}
lon0, lon1, lat0, lat1 = 90, 290, -15, 15
glon = np.arange(lon0, lon1 + 0.01, 1.0); glat = np.arange(lat0, lat1 + 0.01, 1.0); LON, LAT = np.meshgrid(glon, glat)
pts = np.c_[c['lon'], c['lat']]
def grid(v, mask=None):
    m = np.isfinite(v) if mask is None else (mask & np.isfinite(v))
    tri = mtri.Triangulation(pts[m, 0], pts[m, 1])
    return mtri.LinearTriInterpolator(tri, v[m])(LON, LAT).filled(np.nan)
landg = grid(c['land']) > 0.5
polys = d.get('polys')
def cells(ax, values, levels, cmap, mask=None):
    if polys is None: return None
    v = np.array(values, dtype=float); ok = np.isfinite(v) if mask is None else (mask & np.isfinite(v))
    cm = plt.get_cmap(cmap); norm = BoundaryNorm(levels, cm.N, extend='both')
    pc = PolyCollection([polys[i] for i in np.flatnonzero(ok)], array=v[ok], cmap=cm, norm=norm, edgecolors='face', linewidths=0.3, zorder=1)
    ax.add_collection(pc)
    ax.add_collection(PolyCollection([polys[i] for i in np.flatnonzero(~ok)], facecolors='#3a3a3a', edgecolors='#3a3a3a', linewidths=0.3, zorder=2))
    return pc
def smooth(F, passes=6):
    G = F.copy()
    for _ in range(passes):
        P = np.pad(G, 1, mode='edge'); acc = np.zeros_like(G); w = np.zeros_like(G)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                Q = P[1+dy:1+dy+G.shape[0], 1+dx:1+dx+G.shape[1]]; ok = np.isfinite(Q); acc[ok] += Q[ok]; w[ok] += 1
        G = np.where(w > 0, acc / np.maximum(w, 1), np.nan)
    G[~np.isfinite(F)] = np.nan
    return G
title = sys.argv[3] if len(sys.argv) > 3 else default_title(d)
fig, axes = plt.subplots(2, 1, figsize=(16.5, 7.6), facecolor='#0d1117')
xt = [90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 290]; xl = ['90E', '100E', '120E', '140E', '160E', '180', '160W', '140W', '120W', '100W', '80W', '70W']
def style(ax, label):
    ax.set_facecolor('#0d1117'); ax.set_xticks(xt); ax.set_xticklabels(xl, color='w'); ax.tick_params(colors='w'); ax.set_xlim(lon0, lon1); ax.set_ylim(lat0, lat1)
    ax.set_ylabel('Latitude', color='w'); [sp.set_color('w') for sp in ax.spines.values()]; ax.set_title(label, color='w', fontsize=10, loc='left')
ax = axes[0]
T = grid(c['tair']); P = grid(c['mslp']); U = grid(c['u10']); V = grid(c['v10'])
levT = np.arange(21, 31.5, 0.5)
cf = cells(ax, c['tair'], levT, 'RdYlBu_r', mask=c['land'] < 0.5) if polys is not None else ax.contourf(LON, LAT, T, levels=levT, cmap='RdYlBu_r', extend='both')
if polys is None: ax.contourf(LON, LAT, landg.astype(float), levels=[0.5, 1.5], colors=['#3a3a3a'])
cs = ax.contour(LON, LAT, smooth(P, 2), levels=np.arange(990, 1030, 2.5), colors='k', linewidths=0.6, zorder=3); ax.clabel(cs, fmt='%.1f', fontsize=7, inline=True)
WIND_SCALE, CURRENT_SCALE = 220, 9
s = 3; ax.quiver(LON[::s, ::s], LAT[::s, ::s], U[::s, ::s], V[::s, ::s], color='w', edgecolor='k', linewidth=0.5, scale=WIND_SCALE, width=0.0022, headwidth=3.5, headlength=4, zorder=4)
cb = fig.colorbar(cf, ax=ax, pad=0.01, aspect=25); cb.set_label('lowest-layer air temperature °C', color='w'); cb.ax.yaxis.set_tick_params(color='w'); plt.setp(cb.ax.get_yticklabels(), color='w')
style(ax, f"Lowest-layer wind (arrows, 5 m/s ≈ {100 * 5 / WIND_SCALE:.1f}% of the axis width), air temperature (colour) and sea-level pressure (2.5 hPa isobars) — {title}")
ax = axes[1]
sea = c['land'] < 0.5
S = grid(c['sst'], sea); Z = grid(c['thermocline'], sea); CU = grid(c['cu'], sea); CV = grid(c['cv'], sea)
S[landg] = np.nan; Z[landg] = np.nan
levS = np.arange(23, 30.5, 0.5)
cf = cells(ax, c['sst'], levS, 'RdYlBu_r', mask=sea) if polys is not None else ax.contourf(LON, LAT, S, levels=levS, cmap='RdYlBu_r', extend='both')
if polys is None: ax.contourf(LON, LAT, landg.astype(float), levels=[0.5, 1.5], colors=['#3a3a3a'])
cs = ax.contour(LON, LAT, smooth(Z), levels=np.arange(50, 400, 25), colors='k', linewidths=0.7, zorder=3); ax.clabel(cs, fmt='%d m', fontsize=7, inline=True)
mag = np.hypot(CU, CV); ref = 0.3
shrink = np.where(mag > 0, ref * np.sqrt(mag / ref) / np.maximum(mag, 1e-9), 0)
ax.quiver(LON[::s, ::s], LAT[::s, ::s], (CU * shrink)[::s, ::s], (CV * shrink)[::s, ::s], color='w', edgecolor='k', linewidth=0.5, scale=CURRENT_SCALE, width=0.0022, headwidth=3.5, headlength=4, zorder=4)
cb = fig.colorbar(cf, ax=ax, pad=0.01, aspect=25); cb.set_label('mixed-layer temperature °C', color='w'); cb.ax.yaxis.set_tick_params(color='w'); plt.setp(cb.ax.get_yticklabels(), color='w')
style(ax, f"Mixed-layer current (arrows on a square-root length scale: 0.3 m/s ≈ {100 * ref / CURRENT_SCALE:.1f}%, 0.03 m/s ≈ {100 * ref * np.sqrt(0.1) / CURRENT_SCALE:.1f}% of the axis width), mixed-layer temperature (colour) and the top of the {d.get('thermoclineDensity', 1024.0):.1f} class (25 m contours)")
ax.set_xlabel('Longitude', color='w')
plt.tight_layout(); fig.savefig(sys.argv[2], dpi=90, facecolor=fig.get_facecolor()); print('wrote', sys.argv[2])
