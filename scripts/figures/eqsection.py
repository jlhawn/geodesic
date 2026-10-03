# The contoured equatorial Pacific temperature section from eqsection.mjs.
#   python3 scripts/figures/eqsection.py <json> <png> [title]
import sys
import numpy as np
import matplotlib.pyplot as plt
from plotting import load, default_title

d = load(sys.argv[1])
lons = np.array(d['lons']); depths = np.array(d['depths'])
T = np.array([[np.nan if v is None else v for v in col] for col in d['T']]).T
title = sys.argv[3] if len(sys.argv) > 3 else default_title(d)
fig, ax = plt.subplots(figsize=(11, 4.6), facecolor='#0d1117'); ax.set_facecolor('#3a3a3a')
levels = np.arange(10, 31, 1)
cf = ax.contourf(lons, depths, T, levels=levels, cmap='jet', extend='both')
cs = ax.contour(lons, depths, T, levels=levels, colors='k', linewidths=0.6)
ax.clabel(cs, levels=[lv for lv in [14, 16, 18, 20, 22, 24, 26, 28] if lv in cs.levels], fmt='%d', fontsize=9, inline=True)
ax.set_xlim(lons[0], lons[-1]); ax.set_ylim(300, 0)
ticks = [120, 150, 180, 210, 240, 280]
ax.set_xticks(ticks); ax.set_xticklabels(['120E', '150E', '180', '150W', '120W', '80W'], color='w')
ax.tick_params(colors='w'); ax.set_ylabel('Depth (m)', color='w'); ax.set_xlabel('Longitude', color='w')
for sp in ax.spines.values(): sp.set_color('w')
ax.set_title(f"Equatorial Pacific temperature, 2S–2N mean, {d.get('layers', '?')}-layer ocean — {title}", color='w', fontsize=10)
cb = fig.colorbar(cf, ax=ax, pad=0.02); cb.set_label('°C', color='w'); cb.ax.yaxis.set_tick_params(color='w'); plt.setp(cb.ax.get_yticklabels(), color='w')
plt.tight_layout(); fig.savefig(sys.argv[2], dpi=90, facecolor=fig.get_facecolor()); print('wrote', sys.argv[2])
