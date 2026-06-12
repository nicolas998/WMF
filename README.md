___
# wmf
___

WMF (Watershed Modeling Framework) is a module to design to work with hydrographic watersheds and the execution of hydrological models. Additionally, it contains tools for data visualization, analysis of variables, and geomorphological analysis.

## Installation

### Pre-built wheels (no compiler needed)

Pre-compiled wheels for Windows and Linux (64-bit) are built by GitHub
Actions for every tagged release and published to PyPI, so on a supported
platform a plain pip install works without any Fortran compiler:

```bash
pip install wmf
```

Wheels for the development branch can also be downloaded from the
[Actions artifacts](https://github.com/nicolas998/WMF/actions/workflows/wheels.yml)
of the `Build wheels` workflow.

### Building from source

Requirements:

- Python >= 3.9 (tested up to 3.12)
- A Fortran compiler (**gfortran**) and a C compiler
  - Linux: `sudo apt install gfortran` (or your distro equivalent)
  - Windows: install [MSYS2](https://www.msys2.org/) and run
    `pacman -S mingw-w64-x86_64-gcc-fortran`, then put `C:\msys64\mingw64\bin`
    on the `PATH` while installing

The package is built with [meson-python](https://mesonbuild.com/meson-python/)
(the old `numpy.distutils` build no longer works on Python >= 3.12), so a
regular pip install compiles the Fortran extensions automatically:

```bash
python -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install .                   # core install
pip install .[geo,calib]        # + netCDF4/rasterio/gdal + deap
```

Note for Windows: `pip install gdal` has no official wheels; grab one from
[cgohlke/geospatial-wheels](https://github.com/cgohlke/geospatial-wheels/releases)
or use conda.

To import the module:

>- from wmf import wmf

To check the installation, run the smoke test (synthetic watershed, no data
files needed):

```bash
python tests/smoke_test.py
```

## Modules
___

WMF is the result of two Fortran modules and one python script.

**cuencas.f90**

This module has the tools to extract streams and watersheds from a given DIR and DEM maps.  Also, it has multiple functions to perform different tasks based on the topology of the watershed. Some of the tasks are:

- Obtain the area, slope, TI, HAND, Width function, Horton Order.
- Extract the network, identify hillslopes and links. 
- Convert raster to the topology of the watershed.  It can also convert variables to hills or vice-versa.

**modelos.f90**

This module has the TETIS model (Velez, 2001) written from scratch. It also contains several sub-models and has functions to write and read binary data related to the input and output variables of the hydrological model.  The sub-models are:

- HydroFlash: A flash flood 1D hydraulic model to obtain flood plains during execution based on the DEM and the results of the hydrological model.

- SHIA-landslide: An adaptation of the landslide model developed by Aristizabal (2014).

- SED: An adaptation of the sediment production model developed for the CASC2D-SED model.

**wmf.py**

This is the base script that merges cuencas.f90 and modelos.f90.  It has defined several classes such as the **SimuBasin** class that could be considered the heart of WMF.  In this module, we have many functionalities done as an interface to the Fortran modules.

## Dependencies

Installed automatically: numpy, scipy, pandas, matplotlib.

Optional (`pip install .[geo,calib]`):

- netCDF4: save/load basins (`SimuBasin.Save_SimuBasin`)
- gdal/osgeo: read raster and vector maps
- rasterio: basin polygon extraction
- deap: NSGA-II calibration
