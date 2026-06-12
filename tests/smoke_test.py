"""End-to-end smoke test for WMF on Python 3.12.

Builds a synthetic V-shaped valley (a channel along the center column
draining south), writes DEM/DIR GeoTIFFs, and runs the core WMF workflow:
raster reading, stream tracing, basin delineation, geomorphology,
plotting, and SimuBasin netCDF save/load round-trip.

Run from anywhere EXCEPT the repository root (the wmf/ source folder would
shadow the installed package):
    .venv/Scripts/python tests/smoke_test.py   <- works too: we chdir first
"""
import os
import sys
import tempfile

# Avoid the repo root shadowing the installed wmf package
if os.path.isdir('wmf') and os.path.isfile(os.path.join('wmf', 'wmf.py')):
    sys.path = [p for p in sys.path if os.path.abspath(p or '.') != os.path.abspath('.')]

import numpy as np
import matplotlib
matplotlib.use('Agg')
from osgeo import gdal, osr
gdal.UseExceptions()

from wmf import wmf

NCOLS, NROWS = 60, 80
DX = 30.0
XLL, YTOP = 440000.0, 700000.0
EPSG = 32618
NODATA = -9999.0
CHANNEL_COL = 30
OUTLET_ROW = 78


def write_tif(path, arr, dtype):
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(path, NCOLS, NROWS, 1, dtype)
    ds.SetGeoTransform((XLL, DX, 0.0, YTOP, 0.0, -DX))
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(EPSG)
    ds.SetProjection(sr.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.SetNoDataValue(NODATA)
    ds.FlushCache()
    ds = None


def cell_xy(col, row):
    return XLL + (col + 0.5) * DX, YTOP - (row + 0.5) * DX


def main():
    tmp = tempfile.mkdtemp(prefix='wmf_smoke_')
    dem_path = os.path.join(tmp, 'dem.tif')
    dir_path = os.path.join(tmp, 'dir.tif')

    # V-shaped valley: channel along CHANNEL_COL flowing south
    rows = np.arange(NROWS, dtype=float)[:, None]
    cols = np.arange(NCOLS, dtype=float)[None, :]
    dem = 100.0 - 0.5 * rows + 2.0 * np.abs(cols - CHANNEL_COL)

    # DIR in GRASS r.watershed convention (what read_map_raster expects
    # with DIRformat='r.watershed'): 8=E, 4=W, 6=S
    dir_grass = np.full((NROWS, NCOLS), 8, dtype=np.int32)   # west side -> E
    dir_grass[:, CHANNEL_COL + 1:] = 4                       # east side -> W
    dir_grass[:, CHANNEL_COL] = 6                            # channel  -> S

    write_tif(dem_path, dem.astype(np.float32), gdal.GDT_Float32)
    write_tif(dir_path, dir_grass, gdal.GDT_Int32)

    print('--- read_map_raster (DEM, DIR) ---')
    DEM, epsg = wmf.read_map_raster(dem_path, isDEMorDIR=True, dxp=DX, noDataP=NODATA)
    DIR, epsg2 = wmf.read_map_raster(dir_path, isDEMorDIR=True, isDIR=True,
                                     DIRformat='r.watershed', dxp=DX, noDataP=NODATA)
    assert epsg == str(EPSG), f'EPSG mismatch: {epsg}'
    assert DEM.shape == (NCOLS, NROWS), f'DEM shape {DEM.shape}'
    # after reclass: keypad directions E=6 / W=4 / S=2
    assert DIR[CHANNEL_COL - 1, 10] == 6 and DIR[CHANNEL_COL + 1, 10] == 4
    assert DIR[CHANNEL_COL, 10] == 2
    print('DEM/DIR read OK, EPSG =', epsg)

    print('--- Stream trace from channel head ---')
    head_x, head_y = cell_xy(CHANNEL_COL, 5)
    st = wmf.Stream(head_x, head_y, DEM, DIR)
    print('stream cells:', st.structure.shape[1])
    assert st.structure.shape[1] >= OUTLET_ROW - 5, 'stream too short'

    print('--- Basin delineation ---')
    out_x, out_y = cell_xy(CHANNEL_COL, OUTLET_ROW)
    basin = wmf.Basin(out_x, out_y, DEM, DIR, name='synthetic', threshold=300, stream=st)
    expected = NCOLS * (OUTLET_ROW + 1)
    print(f'ncells = {basin.ncells} (expected {expected})')
    assert basin.ncells == expected, 'unexpected basin size'

    print('--- Geomorphology parameters ---')
    basin.GetGeo_Cell_Basics()
    basin.GetGeo_Parameters()
    area = float(basin.GeoParameters.loc['Area_[km2]', 'value'])
    print('Area [km2] =', round(area, 4))
    assert abs(area - expected * DX * DX / 1e6) < 0.05, 'area mismatch'

    print('--- Hypsometric curve + plot (Agg) ---')
    basin.GetGeo_Ppal_Hipsometric(threshold=300)
    png = os.path.join(tmp, 'hipso.png')
    basin.Plot_Hipsometric(path=png, ventana=10)
    assert os.path.getsize(png) > 5000
    print('plot written:', png)

    print('--- SimuBasin + netCDF save/load round-trip ---')
    sim = wmf.SimuBasin(out_x, out_y, DEM, DIR, name='synthetic', threshold=300,
                        dt=300, modelType='cells')
    nc = os.path.join(tmp, 'basin_sim.nc')
    sim.Save_SimuBasin(nc)
    sim2 = wmf.SimuBasin(path=nc)
    print(f'saved ncells = {sim.ncells}, reloaded ncells = {sim2.ncells}')
    assert sim2.ncells == sim.ncells, 'netCDF round-trip mismatch'

    print()
    print('ALL SMOKE TESTS PASSED on Python %d.%d.%d' % sys.version_info[:3])


if __name__ == '__main__':
    main()
