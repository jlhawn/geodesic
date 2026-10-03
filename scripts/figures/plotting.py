# Shared by the figure plots: the dark style, the default title and a raster of cell values on a
# regular longitude–latitude grid, each pixel taking the value of its nearest cell.
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

BACKGROUND = '#0b1020'
LAND = '#3a3a3a'
TEXT = '#e6e9f0'
MUTED = '#9aa3b5'


def load(path):
    with open(path) as f:
        return json.load(f)


def default_title(d):
    return f"{d.get('tag', 'N=' + str(d['N']))} day {d['day']} ({d['season']})" if 'season' in d else f"N={d['N']} day {d['day']}"


def array(values):
    return np.array([np.nan if v is None else v for v in values], dtype=float)


class Raster:
    """Nearest-cell sampling of cell values onto a res-degree grid covering the globe (or a window)."""

    def __init__(self, lon, lat, res, window=(-180, 180, -90, 90)):
        lon = np.asarray(lon, float); lat = np.asarray(lat, float)
        lon = np.where(lon > 180, lon - 360, lon)
        x0, x1, y0, y1 = window
        n = len(lon)
        wrap = np.abs(lon) > 170
        xs = [lon, lon[wrap] - 360 * np.sign(lon[wrap])]; ys = [lat, lat[wrap]]; owner = [np.arange(n), np.nonzero(wrap)[0]]
        caps = np.arange(-190, 191, 2.0)
        for pole in (90, -90):
            nearest = int(np.argmax(lat * np.sign(pole)))
            xs.append(caps); ys.append(np.full(len(caps), float(pole))); owner.append(np.full(len(caps), nearest))
        tri = mtri.Triangulation(np.concatenate(xs), np.concatenate(ys))
        owner = np.concatenate(owner)
        self.gx = np.arange(x0 + res / 2, x1, res); self.gy = np.arange(y0 + res / 2, y1, res)
        LON, LAT = np.meshgrid(self.gx, self.gy)
        t = tri.get_trifinder()(LON, LAT); V = tri.triangles[np.maximum(t, 0)]
        dist = ((tri.x[V] - LON[..., None]) * np.cos(np.radians(LAT))[..., None]) ** 2 + (tri.y[V] - LAT[..., None]) ** 2
        self.cell = owner[np.take_along_axis(V, dist.argmin(-1)[..., None], -1)[..., 0]]
        self.outside = t < 0
        self.extent = (x0, x1, y0, y1)
        self.LON, self.LAT = LON, LAT

    def __call__(self, values):
        G = np.asarray(values, float)[self.cell]
        G[self.outside] = np.nan
        return G


def style_axes(ax, fig=None):
    ax.set_facecolor(BACKGROUND)
    ax.tick_params(colors=MUTED, labelcolor=MUTED, labelsize=8, length=3)
    for sp in ax.spines.values():
        sp.set_color('#444')


def world_ticks(ax):
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-90, 91, 30))
    ax.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    ax.set_yticklabels(['90S', '60S', '30S', 'EQ', '30N', '60N', '90N'])


def colourbar(fig, mappable, ax, label=None, **kwargs):
    cb = fig.colorbar(mappable, ax=ax, **kwargs)
    cb.ax.yaxis.set_tick_params(color=MUTED, labelsize=8)
    plt.setp(cb.ax.get_yticklabels(), color=MUTED)
    cb.outline.set_edgecolor('#444')
    if label:
        cb.set_label(label, color=TEXT, fontsize=9)
    return cb
