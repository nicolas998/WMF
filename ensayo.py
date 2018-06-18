from wmf import wmf 

DEM,mierda = wmf.read_map_raster('/home/nicolas/Dropbox/Trabajos/HZ/05_Buey/raster/demFin.tif', isDEMorDIR=True,
    dxp = 53., noDataP = -9999)
DIR,mierda = wmf.read_map_raster('/home/nicolas/Dropbox/Trabajos/HZ/05_Buey/raster/dirFin.tif', isDEMorDIR=True, 
    isDIR= True, dxp = 53., noDataP = -9999)
    
#st = wmf.Stream(-75.5512, 5.8174, DEM, DIR)

#cu = wmf.SimuBasin(-75.5538, 5.7580, DEM, DIR, stream=st, name='buey', umbral=70)

#print cu.CellSlope
#cu.Save_SimuBasin('Ensayo.nc')
#cu = wmf.SimuBasin(rute='Ensayo.nc')

#print epsg
print wmf.Global_EPSG
