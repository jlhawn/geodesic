# Sources and licenses of the aerosol climatology file

The file is `data/aerosol_mac_2deg.bin` (1,555,624 bytes, SHA-256 `bd707ab14ae9e52432e2dc5447ac6950dd2e61cf44601a13dd2c4ec9344545f7`), written by `scripts/packAerosol.py`. It is the aerosol of every run with `radiation.aerosol: 'mac'`, the default, and the page fetches it when it starts one. The 120 spectral constants of `js/physics/aerosolSpectra.module.js`, written by `scripts/aerosolFit.mjs`, are derived from the same source.

## MACv2

Kinne, S. (2019), The MACv2 aerosol climatology, Tellus B 71, 1623639, doi:10.1080/16000889.2019.1623639. The climatology is Stefan Kinne's, of the Max Planck Institute for Meteorology (MPI-M), Hamburg, for the reference year 2005. Its radiative effects are described in Kinne, S. (2019), Aerosol radiative effects with MACv2, Atmos. Chem. Phys. 19, 10919 to 10959.

The files are MPI-M's, from its anonymous FTP server, ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018/for_ICON/ (the README there is `README_MACv2`). The file is built from four of them, which were fetched on October 26 and 27, 2026:

| file | bytes | modified on the server | SHA-256 |
|---|---|---|---|
| `spec_2005/by_wave/gt_f_00550nm.nc` | 10,112,312 | Aug 13, 2018 | `1b4b958799aa6e47744cc9d9f5f57d59e4543a00007c22300ea08cddf5db01f0` |
| `spec_2005/by_wave/gt_c_00550nm.nc` | 10,112,312 | Aug 13, 2018 | `0eb86e1791656aa2dcd2b19d93a55b4bdf610f8d49ddb7de69327b0e89eb7a70` |
| `altitude/aeropt_fin_alt_km20.nc` | 62,470,752 | Oct 5, 2017 | `4d01453fe3b2b0b57af17e88cf9b3d55b45c5a5ce31ee6a81653f2303cf987f8` |
| `altitude/aeropt_coa_alt_km20.nc` | 62,470,756 | Oct 5, 2017 | `787983ab57ff1f4f15a32de23e3bcc654fe5ef7a6aab959f8c6690c87f6811e2` |

- The two 550 nm files hold the fine and the coarse mode's monthly aerosol optical depth (`aod`), single scattering albedo (`ssa`) and asymmetry (`asy`) on a 1 degree grid; "created by stefan_kinne in 2017_11" is in their history attribute.
- The two altitude files hold, for each month, each mode's share of the column depth in 20 layers above sea level (`zkm_aer_fine_mo` and `zkm_aer_coarse_mo`) and the surface altitude of each cell (`asl`); "created by stefan_kinne in 2016_11" is in their history attribute.
- For the spectral constants, the files of the same names at the other solar wavelengths, `spec_2005/by_wave/gt_f_WWWWWnm.nc` and `gt_c_WWWWWnm.nc` for WWWWW of 00400, 00700, 01000, 01270, 01460, 01780, 02050, 02320, 02790, 03470 and 08000 (10,112,312 bytes each, modified on the server on Aug 13, 2018): their August records.

## What the derivation did

- Resampled the 1 degree grid, 360 by 180 cells, to 2 degrees, 180 by 90: the depth by area, the single scattering albedo (as 1 - SSA) by the depth, the asymmetry by the depth times the SSA, the height by the depth.
- Stored eight fields for each of the twelve months and each cell, one byte each: the 550 nm depth of the fine and of the coarse mode (0 to 1.2, on a square-root scale), their 1 - SSA (0 to 0.3), their asymmetry (0.45 to 0.85) and the mean height of each mode's aerosol above the ground (0 to 8 km). The layout is in the header of `scripts/packAerosol.py`.
- Took the mean height of a mode in a month and a source cell as the mean of the layers' mid-heights weighted by that month's layer shares, less the cell's surface altitude (`asl`), and not below 0.
- Fitted, for the spectral constants, each mode's August depth, 1 - SSA and asymmetry at the other eleven wavelengths to its 550 nm values by least squares over the cells: a ratio for the depth, a line for the other two.

The file is reproducible: `python3 -I scripts/packAerosol.py RAW_DIR data/aerosol_mac_2deg.bin` with the four files above in RAW_DIR gives the same bytes. Every field of every month was also recomputed from the four files with a second implementation (in JavaScript, from the description above), and agrees byte for byte.

## License

MPI-M states no license or terms of use for these files: not on the server, not in `README_MACv2`, not in either paper. The derived file is distributed with the repository and with the page with the attribution above (the project owner's decision of October 27, 2026), and the original data remain MPI-M's and Stefan Kinne's. The derivation, its errors included, is this project's, not MPI-M's, and is not endorsed by MPI-M or the author. If MPI-M states terms for the data, they apply to the file.
