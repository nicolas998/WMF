import os.path

from qgis.core import QgsRasterLayer, QgsMapLayerRegistry

from wmf import wmf

def cargar_mapa_raster (pathMapaRaster):

    retornoCargaLayerMapaRaster = False

    pathMapaRaster = pathMapaRaster.strip ()

    if (os.path.exists (pathMapaRaster)):

        baseNameMapaRaster   = os.path.basename (pathMapaRaster)
        layerMapaRaster = QgsRasterLayer (pathMapaRaster, baseNameMapaRaster)
        QgsMapLayerRegistry.instance ().addMapLayer (layerMapaRaster)

        retornoCargaLayerMapaRaster = layerMapaRaster.isValid ()

    return retornoCargaLayerMapaRaster

def cargar_mapa_dem_wmf (pathMapaDEM, dxp):

    retornoCargaLayerMapaRaster = False

    pathMapaDEM = pathMapaDEM.strip ()

    try:

        demWMF = wmf.read_map_raster (pathMapaDEM, isDEMorDIR = True, dxp = dxp, noDataP = -9999)
        retornoCargaLayerMapaRaster = True

    except:

        retornoCargaLayerMapaRaster = False

    return retornoCargaLayerMapaRaster

