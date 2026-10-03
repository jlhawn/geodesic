# Four global maps of the mixed-layer stratocumulus deck from mlmdeck.mjs, with its box statistics.
#   python3 scripts/figures/mlmdeck.py <json> <png> [title]
import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, PowerNorm
from plotting import default_title
d = json.load(open(sys.argv[1])); cols = d['columns']; A = np.array([[np.nan if v is None else v for v in r] for r in d['rows']], dtype=float)
c = {k: A[:, i] for i, k in enumerate(cols)}
title = sys.argv[3] if len(sys.argv) > 3 else default_title(d)
lon, lat, area, land, ice, on = c['lon'], c['lat'], c['area'], c['land'] > 0.5, c['ice'], c['active'] > 0.5
sea = ~land; open_ = np.where(sea, area * (1 - ice), 0.0)
res = 0.25 if d['N'] > 64 else 0.5
wrap = np.abs(lon) > 170
tri = mtri.Triangulation(np.r_[lon, lon[wrap] - 360 * np.sign(lon[wrap])], np.r_[lat, lat[wrap]]); owner = np.r_[np.arange(len(lon)), np.nonzero(wrap)[0]]
gx = np.arange(-180 + res / 2, 180, res); gy = np.arange(-90 + res / 2, 90, res); LON, LAT = np.meshgrid(gx, gy)
t = tri.get_trifinder()(LON, LAT); V = tri.triangles[np.maximum(t, 0)]
dist = ((tri.x[V] - LON[..., None]) * np.cos(np.radians(LAT))[..., None]) ** 2 + (tri.y[V] - LAT[..., None]) ** 2
cell = owner[np.take_along_axis(V, dist.argmin(-1)[..., None], -1)[..., 0]]
def grid(v):
    G = v[cell].astype(float); G[t < 0] = np.nan
    return G
ss_lon, ss_lat = d['subsolar']
mu = np.sin(np.radians(LAT)) * np.sin(np.radians(ss_lat)) + np.cos(np.radians(LAT)) * np.cos(np.radians(ss_lat)) * np.cos(np.radians(LON - ss_lon))
cat = grid(np.where(land, 2, np.where(ice >= 0.5, 1, 0)))
BOXES = [('SE Pacific', -100, -80, -30, -10), ('California', -135, -120, 20, 35), ('Namibia', 5, 15, -30, -10), ('Canaries', -30, -15, 15, 35), ('Peru', -90, -75, -20, -5),
         ('Australia W', 100, 115, -35, -20), ('N Pacific', 160, 200, 30, 50), ('S Ocean', -180, 180, -60, -45)]
def inbox(b):
    _, x0, x1, y0, y1 = b; x = np.where(lon < x0, lon + 360, lon)
    return (x >= x0) & (x <= x1) & (lat >= y0) & (lat <= y1)
cloud = LinearSegmentedColormap.from_list('cloud', ['#1f4f8f', '#5b95cf', '#b9d6f0', '#f7faff'])
height = LinearSegmentedColormap.from_list('height', ['#6b3410', '#c0692a', '#f0b27a', '#fff2e0'])
fig, axes = plt.subplots(2, 2, figsize=(16, 9.6), facecolor='#0b1020')
def panel(ax, field, cmap, norm, label, extend='max'):
    ax.set_facecolor('#0b1020')
    ax.imshow(cat, extent=(-180, 180, -90, 90), origin='lower', cmap=ListedColormap(['#172133', '#3b4a60', '#3a3a3a']), vmin=-0.5, vmax=2.5, interpolation='nearest')
    im = ax.imshow(np.ma.masked_invalid(field), extent=(-180, 180, -90, 90), origin='lower', cmap=cmap, norm=norm, interpolation='nearest')
    ax.contour(LON, LAT, mu, levels=[0], colors='w', linewidths=0.8, linestyles=':', alpha=0.7)
    ax.plot(ss_lon, ss_lat, marker='+', color='#ffd84d', ms=9, mew=1.6)
    for y in (-60, 60): ax.axhline(y, color='#8892a6', lw=0.5, ls='--', alpha=0.5)
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-90, 91, 30))
    ax.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'], color='#bbb', fontsize=8); ax.set_yticklabels(['90S', '60S', '30S', 'EQ', '30N', '60N', '90N'], color='#bbb', fontsize=8)
    ax.tick_params(colors='#666', length=3); [sp.set_color('#444') for sp in ax.spines.values()]
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, extend=extend); cb.ax.yaxis.set_tick_params(color='#ddd'); plt.setp(cb.ax.get_yticklabels(), color='#ddd', fontsize=8); cb.outline.set_edgecolor('#444')
    ax.set_title(label, color='#eee', fontsize=10.5, loc='left')
thick, cover, h, lwp = (np.where(on, c[k], np.nan) for k in ('thickness', 'cover', 'h', 'lwp'))
act = on & np.isfinite(thick)
panel(axes[0, 0], grid(thick), cloud, PowerNorm(0.5, 0, 600), f"(a) Cloud-layer thickness h − z_b, m, host z_b without sunlight (max {np.nanmax(thick):.0f} m; √ scale)")
for b in BOXES[:-1]:
    _, x0, x1, y0, y1 = b; x0, x1 = (x0 - 360, x1 - 360) if x0 > 180 else (x0, x1)
    for off in ((0,) if x1 <= 180 else (0, -360)): axes[0, 0].add_patch(plt.Rectangle((x0 + off, y0), x1 - x0, y1 - y0, fill=False, ec='#ffd84d', lw=0.7, alpha=0.8))
    axes[0, 0].text(-178 if x1 > 180 else x0, y0 - 6.5 if b[0] == 'SE Pacific' else y1 + 1.5, b[0], color='#ffd84d', fontsize=6.5, alpha=0.9, ha='left')
panel(axes[0, 1], grid(cover), cloud, Normalize(0, 1), "(b) Deck cover: 1 coupled → 0.3 decoupled; radiation uses cover × (1 − sea ice)", extend='neither')
panel(axes[1, 0], grid(h), height, Normalize(0, 2000), f"(c) Inversion height h, m (max {np.nanmax(h):.0f} m)")
panel(axes[1, 1], grid(lwp), cloud, PowerNorm(0.5, 0, 300), f"(d) Liquid water path, g/m² (max {np.nanmax(lwp):.0f}; square-root scale); radiation caps it at 150")
def mean(v, w, m): return np.nansum(v[m] * w[m]) / np.sum(w[m]) if np.sum(w[m]) > 0 else np.nan
def local_time(x0, x1): return ((0.5 * (x0 + x1) - ss_lon) / 15 + 12) % 24
trop = np.abs(lat) <= 60
glob = (f"Ice-free sea with an active deck: {100 * mean(on.astype(float), open_, sea):.1f}% of the globe's, {100 * mean(on.astype(float), open_, sea & trop):.1f}% within 60S–60N; "
        f"deck-covered {100 * mean(np.nan_to_num(c['deckf']) / np.maximum(1e-9, 1 - ice), open_, sea):.1f}%. Area-weighted over active columns: thickness {mean(thick, area, act):.0f} m, "
        f"cover {mean(cover, area, act):.2f}, LWP {mean(lwp, area, act):.0f} g/m², h {mean(h, area, act):.0f} m, cloud base {mean(c['cloudBase'], area, act):.0f} m.")
parts = []
for b in BOXES:
    m = inbox(b) & sea; ma = m & act
    stats = f"{mean(thick, area, ma):.0f} m · {mean(cover, area, ma):.2f} · {mean(lwp, area, ma):.0f} g/m²" if ma.any() else 'no deck'
    parts.append(f"{b[0]} {100 * mean(on.astype(float), open_, m):.0f}% · {stats}" + ('' if b[0] == 'S Ocean' else f" · {local_time(b[1], b[2]):04.1f}h"))
boxes = "Background: land grey · sea ice ≥ 50% slate · sea without a deck navy.   Boxes — active % of ice-free sea · thickness · cover · LWP (area-weighted over active cells) · local solar time:\n  " + ' | '.join(parts[:4]) + '\n  ' + ' | '.join(parts[4:])
print('deck:', glob); print('deck boxes:', ' | '.join(parts))
gate = c['gate']
for b in BOXES:
    m = inbox(b) & sea
    print(f"  {b[0]:12s} gates (% of ice-free sea): " + ', '.join(f"{name} {100 * mean((gate == g).astype(float), open_, m):.0f}" for g, name in enumerate(d['gates']) if g and np.any(m & (gate == g))))
fig.suptitle(f"Marine stratocumulus deck of the mixed-layer model — {title}\nthe deck after four steps from the state, subsolar point (+) {ss_lat:.1f}°, {ss_lon:.1f}°, terminator dotted", color='#eee', fontsize=12)
fig.text(0.01, 0.012, glob + '\n' + boxes, color='#cfd6e4', fontsize=8.3, ha='left', va='bottom', family='DejaVu Sans')
fig.subplots_adjust(left=0.03, right=0.975, top=0.9, bottom=0.135, hspace=0.2, wspace=0.1); fig.savefig(sys.argv[2], dpi=130 if d['N'] > 64 else 90, facecolor=fig.get_facecolor()); print('wrote', sys.argv[2])
