#!/usr/bin/env python3
"""Writes data/woa_annual_1deg.bin, the ocean climatology a fresh start
takes its interior from (js/ocean/climatology.module.js reads it).

Source: the World Ocean Atlas 2023 (NOAA National Centers for Environmental
Information; public domain), annual mean ("00"), all decades ("decav"),
objectively analysed fields (t_an, s_an) on the 1.00 degree grid:

  https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav/1.00/woa23_decav_t00_01.nc
  https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/salinity/netcdf/decav/1.00/woa23_decav_s00_01.nc

Locarnini et al. (2023), WOA23 Volume 1: Temperature, NOAA Atlas NESDIS 89,
https://doi.org/10.25923/54bh-1613; Reagan et al. (2023), WOA23 Volume 2:
Salinity, NOAA Atlas NESDIS 90, https://doi.org/10.25923/70qt-9574.

Usage (numpy and netCDF4 in a virtual environment):

  python3 -m venv venv && venv/bin/pip install numpy netCDF4
  venv/bin/python scripts/packWoa.py woa23_decav_t00_01.nc woa23_decav_s00_01.nc [out.bin]

The file keeps the atlas's 360 x 180 grid at the DEPTHS below, all of them
WOA standard levels. Temperature is potential temperature referenced to
the surface, from the atlas's in-situ ITS-90 temperature by the UNESCO
algorithm (Fofonoff and Millard 1983: Bryden's 1973 adiabatic lapse rate
integrated by fourth-order Runge-Kutta) at the pressure of each depth
(Saunders 1981); salinity is the atlas's practical salinity. A point is
missing wherever either field is.

Layout, little-endian:
  0   char[4]   'WOA1'
  4   uint32    byte offset of the data
  8   uint16    nLon     10  uint16  nLat     12  uint16  nDepth
  14  int16     missing sentinel
  16  float32   lon0, dLon, lat0, dLat   (degrees; centre of column 0 / row 0)
  32  float32   tOffset, tScale, sOffset, sScale
  48  float32   depths[nDepth]  (m)
      uint16    n, then n bytes of ASCII provenance
      zero padding to a multiple of 8
  data: int16 T[nDepth][nLat][nLon], then int16 S[nDepth][nLat][nLon];
        value = offset + scale * stored, rows south to north, columns
        west to east from lon0.
"""
import struct
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

DEPTHS = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000, 4000, 5000]
T_OFFSET, T_SCALE = 15.0, 0.001
S_OFFSET, S_SCALE = 20.0, 0.001
MISSING = -32768


def adiabatic_gradient(s, t, p):
    """Bryden (1973), degC per dbar; t IPTS-68, p dbar."""
    ds = s - 35.0
    return ((((-2.1687e-16 * t + 1.8676e-14) * t - 4.6206e-13) * p
             + ((2.7759e-12 * t - 1.1351e-10) * ds + ((-5.4481e-14 * t + 8.733e-12) * t - 6.7795e-10) * t + 1.8741e-8)) * p
            + (-4.2393e-8 * t + 1.8932e-6) * ds + ((6.6228e-10 * t - 6.836e-8) * t + 8.5258e-6) * t + 3.5803e-5)


def potential_temperature(s, t90, p, reference=0.0):
    """Fofonoff and Millard (1983); ITS-90 in and out."""
    t = t90 * 1.00024
    h = reference - p
    xk = h * adiabatic_gradient(s, t, p)
    t = t + 0.5 * xk
    q = xk
    p = p + 0.5 * h
    xk = h * adiabatic_gradient(s, t, p)
    t = t + 0.29289322 * (xk - q)
    q = 0.58578644 * xk + 0.121320344 * q
    xk = h * adiabatic_gradient(s, t, p)
    t = t + 1.707106781 * (xk - q)
    q = 3.414213562 * xk - 4.121320344 * q
    p = p + 0.5 * h
    xk = h * adiabatic_gradient(s, t, p)
    return (t + (xk - 2.0 * q) / 6.0) / 1.00024


def pressure(depth, lat):
    """Saunders (1981), dbar."""
    c1 = (5.92 + 5.25 * np.sin(np.radians(lat)) ** 2) * 1e-3
    return ((1 - c1) - np.sqrt((1 - c1) ** 2 - 8.84e-6 * depth)) / 4.42e-6


def encode(values, offset, scale):
    stored = np.round((values - offset) / scale)
    if np.nanmin(stored) <= MISSING or np.nanmax(stored) > 32767:
        raise ValueError(f'values {np.nanmin(values)}..{np.nanmax(values)} outside the int16 range at offset {offset}, scale {scale}')
    return np.where(np.isnan(stored), MISSING, stored).astype('<i2')


def main(t_path, s_path, out_path):
    t_file, s_file = Dataset(t_path), Dataset(s_path)
    depth = t_file.variables['depth'][:].data
    if not np.array_equal(depth, s_file.variables['depth'][:].data):
        raise ValueError('temperature and salinity depths differ')
    levels = [int(np.flatnonzero(depth == d)[0]) for d in DEPTHS]
    lat = t_file.variables['lat'][:].data.astype(float)
    lon = t_file.variables['lon'][:].data.astype(float)
    t = np.ma.filled(t_file.variables['t_an'][0, levels, :, :].astype(float), np.nan)
    s = np.ma.filled(s_file.variables['s_an'][0, levels, :, :].astype(float), np.nan)
    missing = np.isnan(t) | np.isnan(s)
    t[missing] = np.nan
    s[missing] = np.nan
    p = pressure(np.array(DEPTHS, float)[:, None, None], lat[None, :, None])
    theta = potential_temperature(s, t, p)
    source = (f'WOA23 decav annual 1.00 degree t_an and s_an (NOAA NCEI, public domain); '
              f'potential temperature (UNESCO 1983) degC, practical salinity').encode('ascii')
    head = struct.pack('<4sIHHHh', b'WOA1', 0, len(lon), len(lat), len(DEPTHS), MISSING)
    head += struct.pack('<8f', lon[0], lon[1] - lon[0], lat[0], lat[1] - lat[0], T_OFFSET, T_SCALE, S_OFFSET, S_SCALE)
    head += struct.pack(f'<{len(DEPTHS)}f', *DEPTHS)
    head += struct.pack('<H', len(source)) + source
    head += b'\0' * (-len(head) % 8)
    head = head[:4] + struct.pack('<I', len(head)) + head[8:]
    with open(out_path, 'wb') as out:
        out.write(head)
        out.write(encode(theta, T_OFFSET, T_SCALE).tobytes())
        out.write(encode(s, S_OFFSET, S_SCALE).tobytes())
    valid = ~missing
    print(f'{out_path}: {Path(out_path).stat().st_size} bytes, {len(lon)} x {len(lat)} x {len(DEPTHS)}; '
          f'theta {np.nanmin(theta):.3f}..{np.nanmax(theta):.3f} degC, S {np.nanmin(s):.3f}..{np.nanmax(s):.3f}; '
          f'ocean points at the surface {valid[0].sum()}, at 5000 m {valid[-1].sum()}')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else str(Path(__file__).resolve().parent.parent / 'data' / 'woa_annual_1deg.bin'))
