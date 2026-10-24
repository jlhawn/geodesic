# Sources and licenses of the n1024 geography files

The files are:
- `data/geography_n1024_N<N>.bin`, written by `scripts/geographyMask.mjs`;
- `data/subgrid_etopo2022_n1024_N<N>.bin`, written by `scripts/subgridTerrain.mjs` with `MASK=n1024`;
- `data/geography1024.bin`, which is generated and not committed.

They are derived from the sources below.

## GSHHG

Global Self-consistent Hierarchical High-resolution Geography, version 2.3.7 (June 15, 2017), full resolution (`gshhs_f.b`), from https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-bin-2.3.7.zip.

Wessel, P., and W. H. F. Smith (1996), A global, self-consistent, hierarchical, high-resolution shoreline database, J. Geophys. Res., 101(B4), 8741–8743.

The files above carry GSHHG's coastline (its land at 60″, `scripts/landWater.mjs`). As derived data they are distributed under GSHHG's license: the GNU Lesser General Public License, version 3 or later (`docs/COPYING.LESSERv3.gshhg`), which supplements the GNU General Public License, version 3 (`docs/COPYINGv3.gshhg`). GSHHG's permission notice, from its LICENSE.TXT:

> Permission to use, copy, modify, and distribute (with no more than a
> reasonable redistribution fee) this data and its documentation for
> any purpose is hereby granted, provided that the above copyright notice
> appear in all copies, that both that copyright notice and this
> permission notice appear in supporting documentation, and that the name
> of GSHHG not be used in advertising or publicity pertaining to
> distribution of the software without specific, written prior permission.
> The GSHHG package may be included in a bundled distribution of software
> for which a fee may be charged. The University of Hawaii (UH) and the
> National Oceanic and Atmospheric Administration (NOAA) make no
> representations about the suitability of this software for any purpose.
> It is provided "as is" without expressed or implied warranty. It is
> provided with no support and without obligation on the part of UH or
> NOAA, to assist in its use, correction, modification, or enhancement.

## Natural Earth

`ne_10m_lakes.geojson`: the lakes, their levels, and the alkaline lakes that are saline or, in Australia, dry pans. Natural Earth, public domain, https://www.naturalearthdata.com/, through https://github.com/nvkelso/natural-earth-vector.

## ETOPO 2022

NOAA National Centers for Environmental Information (2022), ETOPO 2022 15 Arc-Second Global Relief Model, doi:10.25921/fd45-gt74. The 60″ surface and bed are used; the data are public domain.
