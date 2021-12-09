from wmf import wmf 
import numpy as np 
from scipy.spatial import Voronoi, voronoi_plot_2d
from skimage.morphology import dilation, square
import os
import geopandas as geo
import pandas as pd
import osgeo
import copy
from scipy import stats
import pylab as pl
from IPython.display import HTML, display
#try:
import ee
# Trigger the authentication flow.
ee.Authenticate()
# Initialize the library.
ee.Initialize()
#except:
 #   print('Warning: ee not found in the current kernel, try to install it before using\nshp2ee')


def progress(value, max=100):
    return HTML("""
        <progress
            value='{value}'
            max='{max}',
            style='width: 100%'
        >
            {value}
        </progress>
    """.format(value=value, max=max))

def shp2ee(path_shp, type = 'single'):
  '''Converts a shapefile to an ee Feature or collection of Features.
  Parameters:
    - path_shp: path to the shapefile to convert.
    - type:
        - single: if the vector has only one feature.
        - multiple: if the vector has multiple features
  Returns:
    - ee.Feature for type single or ee.FeatureCollection for type multiple'''
  #Read the polygon into geopandas
  shapefile = geo.read_file(path_shp)
  shapefile.to_crs(4326, inplace = True)
  
  #If single, returns just one ee feature
  if type == 'single':
    return ee.Feature(eval(shapefile.iloc[0:1].to_json())['features'][0])
  #If multiple, returns just an ee collection of features
  elif type == 'multiple':
    features = []
    print('Converting shp to ee FeatureCollection...')
    out = display(progress(0, shapefile.shape[0]), 
                  display_id=True)
    for i in range(shapefile.shape[0]):
        geom = shapefile.iloc[i:i+1,:] 
        jsonDict = eval(geom.to_json()) 
        geojsonDict = jsonDict['features'][0] 
        features.append(ee.Feature(geojsonDict)) 
        out.update(progress(i,shapefile.shape[0]))
    print('Done')
    return ee.FeatureCollection(features)

def get_soils_data():
    '''Obtains the soil data from OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02
    and returns it'''
    return ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02')

def get_land_use_data(date1 = '2019-01-01', date2 = '2020-12-01'):
    '''Obtains the land use data from the USDA/NASS/CDL
    Parameters:
        - date1 and date2: initial and end date of the ImageCollection, it grabs the first in the list.
    Returns:
        - ee.Image with the land use'''
    land_use = ee.ImageCollection('USDA/NASS/CDL').filterDate(date1, date2)
    return land_use.first()

class ghost_preprocess():
    
    def __init__(self, watershed, path_dem, seg_threshold = 400, seg_point_distance = 100,thre_sens = 10,
                 focus_map = None, focus_dict = None):
        '''Defines the class to derive the topology required by ghost (.riv and .mesh files)
        Parameters:
            - watershed: A watershed element obtained with wmf.SimuBasin            
            - seg_threshold: minimum distance to divide channels into segments [m].
            - seg_point_distance: the distance of the perpendicular points to the channel segments [m].
        Optional
            - thre_sens: in meters, the sensitivity to transform from channels to segments (not need to edit) 
            - focus_map: load a raster map with categories that represent different levels 
                of detail in the extraction of the stream segments and the mes0.h
                If there is no focus map, all the behavior is set to default category: 0
            - focus_dict: a dictionary with the categories of the focus_map, in each category it has 
                the properties to determine the level of detail to extract the mesh. The dictionary 
                may contain the following properties:
                - seg_threshold (int)
                - seg_point_dstance (int)
                Example: {'0': {'seg_threshold':300},
                          '1': {'seg_threshold':600},
                          '2': {'seg_threshold':800}}
                The properties that are not specified are obtained from the arguments given to the function 
                and are applied for all the region regardless of the focus_map
        Returns:
            - self.wat: a copy of the watershed element
            - self.x and self.y: the X and Y cooridnates of the watershed element
            - self.links and self.links2: element that identifies where channels are'''
        #Copoy the watershed element and defines if it is going to use a focus map and a focus dict to 
        #determine the level of detail of diferent regions in the watershed
        self.wat = watershed
        if focus_map is not None:
            focus, prop, epsg = wmf.read_map_raster(focus_map)
            self.focus_map = self.wat.Transform_Map2Basin(focus, prop)
            self.focus_dict = focus_dict
            self.focus_river = [0]
        else:    
            self.focus_map = focus_map
            self.focust_dict = None
            self.focus_river = None        
        #Define the properties to define the size of the segments and the distance of the mesh to them
        self.threshold = seg_threshold
        self.seg_point_distance = seg_point_distance        
        #Get basic geomorphology 
        self.wat.GetGeo_Cell_Basics()
        #self.wat.GetGeo_StreamOrder()
        self.wat.GetGeo_StreamOrder(threshold=watershed.threshold)
        #Get the coordinates of each element and the links only 
        self.x, self.y = wmf.cu.basin_coordxy(self.wat.structure, self.wat.ncells)
        self.links = self.wat.CellCauce * self.wat.hills_own
        channels = np.zeros(watershed.ncells)
        channels[self.wat.CellAcum>self.wat.threshold-thre_sens]=1
        self.links2 = channels * watershed.hills_own
        #Read the DEM map of the region 
        self.DEM, self.dem_prop, self.epsg = wmf.read_map_raster(path_dem)
    
    def get_segments_topology(self, correct_downstream_elev = True, epsilon = 0.1,
                              default_min_order = 1):
        '''Using channel2segment obtains the topological connection 
        of the segments.
        Parameters:             
            - correct_downstream_elev: Defines if the code will correct or not the elevation of the 
                segement downstream.
            - epsilon: the difference between segments with same elevation [m].
            - default_min_order (optional): Use this value to determine the minimum horton 
                order of the channels to convert to segments when using focus. 
        Returns:
            - self.river_topology: list of list describing the connections among segments.
            - self.river_centers: list with the centers size: prop -1
            - self.river_length: list with the length size: prop -1'''
        #Initializes the focur_river property
        self.focus_river = []
        #Code starts
        start = 0
        prop = [[0,-999,self.x[-1],self.y[-1],-999,-999,self.wat.CellHorton_Hill[-1],self.wat.ncells]]
        new_dest = [0]
        #self.focus_river = [0]
        print('Defining river segments...')
        out = display(progress(0, self.wat.nhills), 
                  display_id=True)
        for c, dest in enumerate(self.wat.hills[1][::-1]):
            d = new_dest[dest]
            p, x1, y2, start, focus_categories = self.__channel2segments__(c+1, start, d)
            prop.extend(p)
            new_dest.append(start)
            out.update(progress(c, self.wat.nhills))            
            self.focus_river.extend(focus_categories)
        print('Done')
        self.river_topology = prop
        self.__get_segments_center_length__()
        self.__get_segment_sinuosity__()
        if correct_downstream_elev:
            corrected = self.__correct_downstream_elevation__(epsilon)
            return corrected
        # If there is focus using the horton order
        focus_by_order = False
        for k in self.focus_dict.keys():  
            if 'min_order' in self.focus_dict[k]:
                focus_by_order = True
        if focus_by_order:
            self.select_segments_using_focus(default_min_order)
            
    
    def select_segments_using_focus(self, default_min_order = 1):
        '''Function to obtain the segments that are going to be used based in the focus regions'''            
        #Finds the segments that are above the horton order for that focus level
        a = copy.copy(self.river_topology)
        new_seg = [a[0]]
        new_foc = []
        cont = 0
        for seg, foc in zip(a[1:], self.focus_river):
            try:
                min_order = self.focus_dict[str(foc)]['min_order']
            except:
                min_order = default_min_order
            if seg[-2] >= min_order:
                new_seg.append(seg)
                new_foc.append(foc)
            cont+=1

        #Finds the segments that are orphan after erasing based in the focus areas
        dest = np.array(new_seg).T[1].astype(int).tolist()
        ids = np.array(new_seg).T[0].astype(int).tolist()
        bads = []
        new_seg2 = [new_seg[0]]
        self.focus_river = [new_foc[0]]
        cont = 1
        for i,j in zip(dest[1:],ids[1:]):
            if i not in ids or i in bads:
                bads.append(j)
            else:
                new_seg2.append(new_seg[cont])
                self.focus_river.append(new_foc[cont-1])
            cont+=1

        #Reasign the IDs of the segments now that some of them were wiped out
        dest = np.array(new_seg2).T[1].astype(int)
        for cont, seg in enumerate(new_seg2):
            if seg[0] > cont:
                for d in np.where(dest == seg[0])[0].tolist():
                    new_seg2[d][1] = cont
                seg[0] = cont
        self.river_topology = new_seg2

        #Compute again que properties of the segmentes 
        self.__get_segments_center_length__()
        self.__get_segment_sinuosity__()
    
    def __correct_downstream_elevation__(self, epsilon):
        '''Checks if the elevation of the segment downstream is higher or equal,
        if yes, it corrects it by decreasing it by an epsilon'''
        corrected = []
        dest = np.array(self.river_topology).T[1].astype(int)
        elev = np.array(self.river_topology).T[5]
        for fid in range(len(self.river_topology)-1,0,-1):
            pos = np.where(dest == fid)[0]
            _e = elev[pos].tolist()
            _e.append(self.river_topology[fid][5])
            if np.min(_e) < self.river_topology[fid][5]:
                self.river_topology[fid][5] = np.min(_e) - epsilon
                self.river_topology[fid][4] = self.river_topology[fid][5] - 20
                elev = np.array(self.river_topology).T[5]
                corrected.append(fid)
                #print(fid, np.min(_e), self.river_topology[fid][5])
        return corrected
    
    def get_mesh_river_points(self, clean_close_points = True, 
                              min_river2river_distance = 50):
        '''Obtains the two centers of the lines that are perpendicular 
        to each river segment
        Parameters:
            - clean_close_points: if true, remove lower order points that are at a distance 
            less than min_river2river_distance, sho.
            -min_river2river_distance: The minimum distance between two points extracted from the segments [m].'''
        #Determines if it is going to use the focus dicst for this procedure
        use_focus = False
        if len(self.focus_river) == len(self.river_topology):
            use_focus = True
        
        #Obtains the points adjacent to each segment
        Xp = []
        Yp = []
        for l in range(1,len(self.river_topology)):            
            #Get the coordinates
            xo = self.river_topology[l][2]
            yo = self.river_topology[l][3]
            d = self.river_topology[l][1]
            xd = self.river_topology[d][2]
            yd = self.river_topology[d][3]
            #Get the centroid coordinate
            xc, yc = self.river_center[l-1]
            #Compute the slope
            if xo==xd:
                m = 0.0
            elif yo==yc:
                m = -90
            else:
                m = -1 / ((yo-yd)/(xo-xd))
            #Determine the distance of the points 
            if use_focus is True:
                focus_level = str(self.focus_river[l])
                try:
                    seg_point_distance = self.focus_dict[focus_level]['seg_point_distance']
                except:
                    print('Warning: seg_point_distance not set in the self.focus_dict, will use ')
                    print('the program will use the self.seg_point_distance value instead')
                    seg_point_distance = self.seg_point_distance
            else:
                seg_point_distance = self.seg_point_distance
            #Compute the new X and Y
            x1 = xc - seg_point_distance/(np.sqrt(1+m**2))
            y1 = m*(x1-xc)+yc

            x2 = xc + seg_point_distance/(np.sqrt(1+m**2))
            y2 = m*(x2-xc)+yc

            Xp.append([x1,x2])
            Yp.append([y1,y2])
        Xv = np.array(Xp).T.reshape(len(Xp)*2)
        Yv = np.array(Yp).T.reshape(len(Yp)*2)
        XYr = np.vstack([Xv, Yv])
        self.mesh_points_river = XYr
        
        #Clean close points
        if clean_close_points:
            #Determine the dist to erase close points
            if min_river2river_distance > self.seg_point_distance:
                print('Warning: river to river point distance is greater than segement points distance\
                    the program will set it equal to dist*0.5')
                min_river2river_distance = self.seg_point_distance*0.5
            #Clean the close points
            self.__clean_river_points__(min_river2river_distance)
    
    def get_mesh_grid_points(self, mesh_spaces = 8, border_iter = 3,clean_with_river = True, 
                            min_dem2river_distance = 100):
        '''Obtains the X and Y coordinates for the mesh and the borders of the mesh.
        Parameters:
            - mesh_spaces: step to take values in the grid given by the DEM [pixels].
            - border_iter: number of times to perform the dilation that generates the 
                border points (Must be tested).
            - clean_with_river: if True, it will remove mesh points that are below the 
                min_dem2river_distance value.
            - min_dem2river_distance: The minimum distance between grid points extracted from the DEM 
                and the points extracted from the segements with self.get_mesh_river_points [m].
        Returns:
            - self.mesh_points_dem'''
        # Get the X and Y values fom the DEM mesh
        x_vect, y_vect = wmf.cu.basin_coordxy(self.wat.structure, self.wat.ncells)
        x_mask = wmf.cu.basin_float_var2map(self.wat.structure, x_vect, wmf.cu.ncols,wmf.cu.nrows, self.wat.ncells)
        y_mask = wmf.cu.basin_float_var2map(self.wat.structure, y_vect, wmf.cu.ncols,wmf.cu.nrows, self.wat.ncells)

        #Apply the focus algorithm if active
        if self.focus_map is not None:
            print('Extracting points for the mesh by focus areas...')
            #Obtains a map oif the focus areas following the size of the DEM
            focus_map = wmf.cu.basin_float_var2map(self.wat.structure, 
                                                   self.focus_map, wmf.cu.ncols,wmf.cu.nrows, self.wat.ncells)
            #Iterate through the focus categories
            XYm = []
            for k in self.focus_dict.keys():
                try:
                    mesh_spaces_temp = self.focus_dict[k]['mesh_spaces']
                except:
                    print('Warning: %s category does not have a defined mesh_spaces' % k)
                    print(' will use the default value')
                    mesh_spaces_temp = mesh_spaces                    
                x_steps = np.arange(0,x_mask.shape[0],mesh_spaces_temp)
                y_steps = np.arange(0,x_mask.shape[1],mesh_spaces_temp)
                xv,yv = np.meshgrid(x_steps, y_steps)
                for i in x_steps:
                    for j in y_steps:
                        if x_mask[i,j] > 0 and focus_map[i,j] == int(k):
                            XYm.append([x_mask[i,j],y_mask[i,j]])                            
            XYm = np.array(XYm).T
            print('Done')
        else:
            #Obtain a set of points spaced in the grid
            x_steps = np.arange(0,x_mask.shape[0],mesh_spaces)
            y_steps = np.arange(0,x_mask.shape[1],mesh_spaces)
            xv,yv = np.meshgrid(x_steps, y_steps)
            # Get the xy points inside the 
            print('Extracting points for the mesh without focus areas...')
            XYm = []
            for i in x_steps:
                for j in y_steps:
                    if x_mask[i,j] > 0:
                        XYm.append([x_mask[i,j],y_mask[i,j]])
            XYm = np.array(XYm).T
            print('Done')

        print('Creating border elements...')
        borders = []
        for j in [1,5]:
            #Increase the border n times
            new = np.copy(x_mask)
            new[new > 0] = 1    
            for i in range(border_iter*j):
                old = np.copy(new)
                new = dilation(old, square(3))
            border = new - old
            border = np.where(border>0)

            #Converts the border to a list of X and Y points 
            x_border = []
            y_border = []
            for col in border[0]:
                x_border.append(wmf.cu.dxp*(col - 0.5) + wmf.cu.xll)
            for row in border[1]:
                y_border.append(wmf.cu.dxp*(wmf.cu.nrows - row + 0.5) + wmf.cu.yll)
            borders.append(np.vstack([x_border, y_border]))
        print('Done')
        self.mesh_points_dem = XYm 
        self.mesh_points_boundary = borders
        if clean_with_river:
            print('Cleaning mesh points with the points of the river network...')
            self.__clean_mesh_points__(min_dem2river_distance)
            print('Done')
    
    def get_voronoi_polygons(self):
        # Get the array with all the points for the voronoi
        XYall = np.hstack([self.mesh_points_river,self.mesh_points_dem,self.mesh_points_boundary[0]]).T
        categories = np.hstack([np.ones(self.mesh_points_river.shape[1]),
                                np.ones(self.mesh_points_dem.shape[1])*2,
                               np.ones(self.mesh_points_boundary[0].shape[1])*3])
        self.mesh_points_all = XYall
        self.vor = Voronoi(XYall)
        self.vor_cat = categories
    
    def define_polygons_topology(self, define_left_right = True, reduce_exterior_faces = False,
                               n_exterior_faces = 14, min_exterior_distance = 50,
                               zmethod='mean',zpercentile=50, zkernel=5):
        '''obtains self.polygons_topology and self.polygons_topology_cat after self.get_voronoi_polygons()'''
        poly_prop = []        
        #Estimates the number of accepted polygons
        self.polygons_expected_number = self.vor_cat[self.vor_cat < 3].shape[0]+2 
        #Progress bar
        out = display(progress(0, self.vor.points.shape[0]), 
                      display_id=True)
        #Iterate through the polygons to see if they are valid
        cont = 1
        self.polygons_id = []
        self.all_polygons_id = []
        for poly in range(self.vor.points.shape[0]):
            #First check if not from the buffer
            if self.vor_cat[poly] < 3:                
                # Get the properties of the polygon
                _p_prop = list(self.get_polygon_prop(poly,False,zmethod=zmethod, zpercentile=zpercentile,zkernel=zkernel))
                # determine the number of neighbors that are inside the watershed
                good_neighbors = [i for i in _p_prop[4] if i < self.polygons_expected_number]
                # The polygon is valid if it has at least one neighbor that is not part of the buffer
                if len(good_neighbors) > 1:
                    poly_prop.append(_p_prop)
                    self.polygons_id.append(cont)
            self.all_polygons_id.append(cont)
            cont+=1
            #update the progress bar
            out.update(progress(cont-2, self.vor.points.shape[0]))
        #Define and empty list to store the coorditantes of the polygons
        xyp = []
        #Set the total number of good polygons to be stored
        self.polygons_final_number = len(poly_prop)
        #Iterate trhoug the polygons to get the coordinates
        for count, p in enumerate(poly_prop):
            #Get the coordinates of the polygon
            xyp.append(p[0])
            #If the polygon has negative elevation, computes it as the mean of its neighbors
            if p[1] < 0:
                poly_prop[count][1] = np.mean([poly_prop[i][1] for i in p[4] if i < self.polygons_expected_number and poly_prop[i][1]>0])
        #Set the coordinates and the polygons prop as class attributes
        self.polygons_topology = poly_prop 
        self.polygons_xy = np.array(xyp)
        #convert the neighbors to an array 
        self.__neighbors2array__()
        self.__correct_numbering__()
        # Defines the left and right of the river topo
        if define_left_right:
            self.__get_left_right__(self.polygons_xy, self.river_center)
        # Search for polygons with elevation equal to nan and corrects them 
        zmean = [polygon[1] for polygon in self.polygons_topology if polygon[1]>1]
        zmean = np.mean(zmean) # obtains the mean elevation
        # Iterate to correct polygons with wrong elevation
        for c, polygon in enumerate(self.polygons_topology):
            if np.isnan(polygon[1]):
                #Get the elevation of the neighbors
                z = []
                for p in polygon[4]:
                    try:
                        if self.polygons_topology[p-1][1]>0:
                            z.append(self.polygons_topology[p-1][1])
                    except:
                        pass
                #If the elevation of the neighbors is not nan, set the elevation of the polygon to the mean of the neighbors
                if len(z)>0:
                    self.polygons_topology[c][1] = np.mean(z)
                #Otherwise, set the elevation to the mean of the watershed
                else:
                    self.polygons_topology[c][1] = zmean
        #Reduces the number of exterior faces if the user wants to 
        if reduce_exterior_faces:
          self.decrease_num_faces(n_exterior_faces, min_exterior_distance)

    def check_if_neighbors_exist(self):
        '''Checks if any polygon has a neighbor with an ID higher than the total number of elements'''
        #Set the number of elements and an empty list of bad polygons
        npolygons = len(self.polygons_topology)
        bad_neighbors = []
        #Iterate trhough all the polygons
        for c, polygon in enumerate(self.polygons_topology):
            neighbors = np.array(polygon[4])
            existing = np.array(polygon[7])
            neighbors = neighbors * existing
            if neighbors.max() > npolygons:
                bad_neighbors.append([c+1, neighbors.argmax()])
        #Returns a list that has the polygon followed by the polygon Id that is outside of the range
        return bad_neighbors
    
    def decrease_num_faces(self, step = 2, minfaces = 18):
        '''Reduces the number of faces of polygons next to border polygons. The step indicates 
        the step taken to reduce faces (border[::step]), minfaces indicates the minimum number of faces'''
        def compute_distance(x,y):
            d = np.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2)
            return number2dec(d)
        def number2dec(x):
            t = '%.3f' % x
            return float(t)
        def multiple_distances(coor):
            distances = []
            for i,j in zip(coor[:-1], coor[1:]):
                distances.append(compute_distance(i,j))
            distances.append(number2dec(compute_distance(coor[-1],coor[0])))
            return distances

        def compute_area(coord):
            x = []; y = []
            for i in coord:
                x.append(i[0])
                y.append(i[1])
            pgon = Polygon(zip(x, y)) # Assuming the OP's x,y coordinates
            return pgon.area
        
        cont = 0
        #Iterates through all the polygons to check
        for p,n,i in zip(self.polygons_topology ,self.neighbors, self.polygons_id):
            if n.size > 20 and n.min() == 0:                
                nfaces = n.size
                print('polygon %d has %d faces' % (i, nfaces))
                a = p        
                n = n.tolist() #Neighbors ids
                l = [number2dec(i) for i in a[5]] # Lenghts of the faces
                d = a[6] #Distance between polygons
                b = a[7] # Border neighbors
                c = a[-1] # Coordinates of the polygon

                # Find the order of the coordinates that matches the neighbors order
                ltemp = multiple_distances(c) # Compute the distances in the order of the coordinates 
                c2 = []
                for i in ltemp:
                    c2.append(c[l.index(i)])

                #Find where neighbors start to become border
                b2 = np.array(b)
                diff = b2[:-1] - b2[1:]
                pos = diff.tolist().index(1) +1   

                while nfaces < minfaces:
                    nn = n[:pos] + n[pos::step]
                    nfaces = len(nn)
                    if nfaces > minfaces:
                        step += 1    
                nn = n[:pos] + n[pos::step]
                nfaces = len(nn)
                print('N faces reduced to %d using an step of %d' % (nfaces, step))
                #New polygon neighbors, borders, lenght, and etc...
                nn = n[:pos] + n[pos::step]
                ln = l[:pos] + l[pos::step]
                bn = b[:pos] + b[pos::step]
                cn = c2[:pos] + c2[pos::step]
                dn = d[:pos] + d[pos::step]
                area = compute_area(cn)
                self.polygons_topology[cont][2] = area
                self.polygons_topology[cont][3] = len(nn)
                self.polygons_topology[cont][4] = nn
                self.polygons_topology[cont][5] = ln
                self.polygons_topology[cont][6] = dn  
                self.polygons_topology[cont][7] = bn
                self.neighbors[cont] = np.array(nn)
            cont+=1
    
    def get_polygon_prop(self, elem, plot = False,zmethod = 'mean', zpercentile=50, zkernel=5):    
        '''Get the properties of a polygon, 
        The function is called by the self.define_polygons_topology function'''
        # Get the region from the vor element computed by the Delauney triangulation
        region = self.vor.regions[self.vor.point_region[elem]]
        # Extract the coordinates of the polygon vertices
        polygon = [self.vor.vertices[i] for i in region]
        #Extract the centroid of the polygon
        cent = self.vor.points[elem]        
        #Get the elevation of the polygon (may be improved in the near future to use a kernel and the minimum value in it)
        Z = self.__get_z_from_dem__(cent[0],cent[1],kernel=zkernel, method=zmethod, percentile=zpercentile)
        #compute the projected area of the polygon 
        area = self.__polygon_area__(np.array(polygon).T[0],np.array(polygon).T[1])
        #plot the polygon if the user wants to
        if plot:
            if self.vor_cat[elem] < 3:
                pl.fill(*zip(*polygon), zorder = 1, alpha = 0.5)
            else:
                pl.fill(*zip(*polygon), zorder = 1, alpha = 0.5, c = 'k')
            pl.scatter(cent[0],cent[1], c = 'k', zorder = 0)
        #conver the vertices of the region to a set
        list1_as_set = set(region)
        #Get the number of points in the region
        n_points = self.vor.point_region.shape[0]
        #Empty list to store the neighbors of the polygon
        neighbors = []
        is_border = []
        lface = []
        dneigh = []
        #Iterate over all the elements to determine which are neighbors
        for r in range(n_points):
            #Get a region
            region2 = self.vor.regions[self.vor.point_region[r]]
            #See if both regions intersect
            intersection = list1_as_set.intersection(region2)
            intersection_as_list = list(intersection)
            # Eval if the intersection is not empty
            if len(intersection_as_list) > 1 and len(intersection_as_list)<3: #Not sure why minus three
                if self.vor_cat[r] < 3:
                    #Determine if the neighbor is a border
                    is_border.append(1)
                else:
                    #Or not
                    is_border.append(0)
                #Append the region to the list of neighbors
                neighbors.append(r+1)
                #Get the vertices coordinates of the neighbor
                v1 = np.array(self.vor.vertices[intersection_as_list][0])
                v2 = np.array(self.vor.vertices[intersection_as_list][1])
                #Get the centroid of the neighbor
                c1 = self.vor.points[r]
                #Get the face and centroid to centroid distance
                lface.append(np.linalg.norm(v1- v2, ord = 2, axis=0))
                dneigh.append(np.linalg.norm(c1- cent, ord = 2, axis=0))
                #If users wants to plot the polygons
                if plot:
                    pl.scatter(c1[0], c1[1],c = 'r')
        #If user wants to plot the polygons
        if plot:   
            for i in neighbors:
                region2 = self.vor.regions[self.vor.point_region[i]]
                polygon2 = [self.vor.vertices[i] for i in region2]
                if self.vor_cat[i] < 3:
                    pl.fill(*zip(*polygon2), c = 'g', alpha = 0.5)
                else:
                    pl.fill(*zip(*polygon2), c = 'r', alpha = 0.5)
        #returns the centroid, elevation, area, number of faces, list of neighbors, list of face lenghts, list of centroid to centroid distances
        return cent, Z, area, len(lface), neighbors,lface, dneigh, is_border, polygon
    
    def get_physical_prop(self, ee_data, scale = 30, prop_name = 'feature', band = 'b0',
        sliced = False, xdivisions = None, ydivisions = None):
        '''Extract a property from an ee Image and assign it to the mesh polygons of the object.
        Parameters:
            - ee_data: an ee Image with the desired property.
            - band: name of the band in the ee Image to use.
            - prop_name: name of the property to assign to the self.polygon_shp
            - scale: scale of resolution of the image in m.
            - sliced: work the watershed as a whole or split it.
            - xdivisions: Number of divisions to make over the vector to obtain the 
                properties, same for ydivisions
        Returns:
            - updates the self.polygon_shp with the property.'''
        #Define functions for several tiles 
        def get_boundaries(vector):
            xvalues = []
            yvalues = []
            for geometry in vector['geometry']:
                xvalues.append(geometry.bounds[0])
                xvalues.append(geometry.bounds[2])
                yvalues.append(geometry.bounds[1])
                yvalues.append(geometry.bounds[3])
            return (np.min(xvalues), np.max(xvalues),np.min(yvalues), np.max(yvalues))

        def define_slices(boundaries, xdivisions, ydivisions):
            xsteps = np.linspace(boundaries[0], boundaries[1],xdivisions)
            ysteps = np.linspace(boundaries[2], boundaries[3],ydivisions)
            return xsteps, ysteps

        def geopandas2ee(geoData):
            features = []
            for i in range(geoData.shape[0]):
                geom = geoData.iloc[i:i+1,:] 
                jsonDict = eval(geom.to_json()) 
                geojsonDict = jsonDict['features'][0] 
                features.append(ee.Feature(geojsonDict)) 
            return ee.FeatureCollection(features)
        
        #Get the soils information for each polygon 
        if sliced is False:
            d = ee_data.reduceRegions(self.polygon_ee, reducer = ee.Reducer.mode(), scale = scale).getInfo()
            value = []
            for i in range(self.polygons_shp.shape[0]):
                value.append(d['features'][i]['properties'][band])
            self.polygons_shp[prop_name] = value
        else:
            print('Will extract EE features by tiles')
            shp_copy = self.polygons_shp.copy()
            shp_copy.to_crs(4326, inplace = True)
            self.polygons_shp[prop_name] = 0
            xslice, yslice = define_slices(get_boundaries(shp_copy), xdivisions, ydivisions)
            for x1,x2 in zip(xslice[:-1], xslice[1:]):
                for y1,y2 in zip(yslice[:-1], yslice[1:]):
                    #Obtains an slice of the vector
                    sliced = shp_copy.cx[x1:x2,y1:y2]
                    ee_feature = geopandas2ee(sliced)
                    #Get the properties for that slice
                    d = ee_data.reduceRegions(ee_feature, reducer = ee.Reducer.mode(), scale = 30).getInfo()
                    value = []
                    for i in range(sliced.shape[0]):
                        value.append(d['features'][i]['properties'][band])
                    sliced[prop_name] = value
                    print('slice x: %.1f - %.1f and y: %.1f - %.1f done' % (x1,x2,y1,y2))
                #Merge the properties with the global vector
                self.polygons_shp.loc[sliced.index, prop_name] = sliced[prop_name]
    
    def write_mesh_file(self, path, shp_path = None):
        '''Writes the .mesh file for ghost, also writes the .shp file if it is given.'''
        #Opens the file to write the mesh
        f = open(path,'w')
        #Gets the number of polygons
        n_points = len(self.polygons_topology)
        #Print the progress bar
        print('Writing mesh file...')
        out = display(progress(0, n_points), 
                  display_id=True)
        #Write the header
        f.write('NUMELE\t%d\n' % n_points)
        f.write('INDEX   X   Y   Zmin    Zmax    Area    nFaces\n')
        #Wrtites the properties for each polygon
        for c, p in zip(self.polygons_id,self.polygons_topology):
            f.write('%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n' % (c,p[0][0],p[0][1],p[1]-20,p[1],p[2],p[3]))
        #Writes the header information regarding the neighbors
        f.write('\n')
        f.write('INDEX   ID_Neighbors    Lenght_Face_Neighbors   Distance_to_Neighbors\n')
        #Writes the neighbors for each polygon
        for c, p in enumerate(self.polygons_topology):
            f.write('%d ' % (c+1))
            neighbors = self.neighbors[c] #(np.array(p[4])+1 - self.bad_polygons)*np.array(p[-2])
            for n in neighbors:
                f.write('%d ' % n)
            for lf in p[5]:
                f.write('%.2f ' % lf)
            for lf in p[6]:
                f.write('%.2f ' % lf)
            f.write('\n') 
            out.update(progress(c, n_points))
        f.close()
        #End of mesh file
        print('Mesh file written')
        #If given, write the shp file
        if shp_path is not None:
            print('writing shapefile...')
            #Function to write the shp file
            self.__write_mesh_shp__(shp_path)
            print('Done')
            #Try to generate the ee version of the shapefile to use it in google earth engine
            self.polygons_shp = geo.read_file(shp_path)            
            try:
                self.polygon_ee = shp2ee(shp_path, type='multiple')
            except:
                print('Warning: self.polygon_ee not defined it seems that you dont have ee set up.')
    
    def write_attribute_file(self, path, name_1, name_2):
        '''Writes the att file for ghost using the self.polygons_shp variable and 
        its columns describing the soils texture and land use.
        Parameters:
            - path: out path to write the att file.
            - name_1: name of the column with the soils.
            - name_2: name of the column with the land use.
        Returns:
            - Null, it just writes the att.'''
        f = open(path, 'w')
        f.write('INDEX SOIL LC METEO LAI SS LAKE CLOSE_SEG BCns\n')
        for i in self.polygons_shp.index:
            poly = self.polygons_shp.loc[i,'polygon']
            soil = self.polygons_shp.loc[i,name_1]
            land = self.polygons_shp.loc[i,name_2]
            sides = self.polygons_topology[i][3]
            f.write('%d %d %d 1 1 0 0 1 ' % (poly,soil, land))
            for z in range(sides):
                    f.write('0 ')
            f.write('\n')
            #self.attr !this is for sam
        f.close()
    
    def write_river_file(self, path, shp_path = None):
        prop = self.river_topology
        f = open(path, 'w')
        f.write('NUMRIV\t%d\n' % len(prop[1:]))
        f.write('INDEX	X	Y	ZMIN	ZMAX	LENGTH	DOWN	LEFT	RIGHT	SHAPE	MATRL	BC	RES	XAREA	INACT	LAKE	LRIV\n')
        prop[1][1] = -4
        for seg,cen,lr,lriv,len_riv in zip(prop[1:],self.river_center,self.river_left_right,self.river_sinuosity, self.river_length):
            f.write('%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t' % (seg[0], cen[0], cen[1], seg[4],seg[5],len_riv,seg[1]))
            f.write('%d\t%d\t%d\t%d\t' % (lr[0],lr[1], seg[-2], seg[-2]))
            for i in range(5):
                f.write('0\t')
            f.write('%.2f\n' % lriv)
        max_h = int(np.array(self.river_topology).T[-2].max())
        f.write('SHAPE\t%d\n' % max_h)
        f.write('INDEX\tDPTH\tOINT\tCWID\n')
        for i in range(1,max_h+1):
            width = i**3.3
            width = '%d' % width
            f.write('%d\t0.0 1 %s\n' % (i,width))
        f.write('MATERIAL\t%d\n' % max_h)
        f.write('INDEX\tROUGH\tCWR\tHWEIR\tKEFF\tHSED\n')
        for i in range(1,max_h+1):
            f.write('%d\t0.06\t0.010\t0.20\t2.00E-09\t0.1\n' % i)
        f.write('BC\t0\nRES\t0')
        f.close()
        if shp_path is not None:
            self.__write_river_shp__(shp_path)
            
    def __write_river_shp__(self, path):
        prop = self.river_topology
        DriverFormat='ESRI Shapefile'    
        sr = osgeo.osr.SpatialReference()
        sr.ImportFromEPSG(int(self.wat.epsg))
        driver = osgeo.ogr.GetDriverByName(DriverFormat)
        if os.path.exists(path):
            driver.DeleteDataSource(path)
        shapeData = driver.CreateDataSource(path)
        layer = shapeData.CreateLayer('layer1', sr, osgeo.ogr.wkbLineString)
        layerDefinition = layer.GetLayerDefn()
        names = ['segment','long[m]','z[m]','down','left','right','horder']
        fmts = [osgeo.ogr.OFTInteger,osgeo.ogr.OFTReal,osgeo.ogr.OFTReal,
               osgeo.ogr.OFTInteger, osgeo.ogr.OFTInteger, osgeo.ogr.OFTInteger,
               osgeo.ogr.OFTInteger]
        for name,fmt in zip(names, fmts):
            new_field=osgeo.ogr.FieldDefn(name,fmt)
            layer.CreateField(new_field)
        featureFID = 0
        prop[1][1] = 0
        for l in range(1,len(self.river_topology)):
            xo = prop[l][2]
            yo = prop[l][3]
            d = prop[l][1]
            xd = prop[d][2]
            yd = prop[d][3]
            line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
            line.AddPoint_2D(float(xo),float(yo))
            line.AddPoint_2D(float(xd),float(yd))
            feature = osgeo.ogr.Feature(layerDefinition)
            feature.SetGeometry(line)
            feature.SetFID(0)
            feature.SetField('segment',int(prop[l][0]))
            feature.SetField('long[m]',float(self.river_length[l-1]))
            feature.SetField('z[m]',float(prop[l-1][5]))
            feature.SetField('down',int(d))
            feature.SetField('left',int(self.river_left_right[l-1][0]))
            feature.SetField('right',int(self.river_left_right[l-1][1]))
            feature.SetField('horder',int(prop[l][-2]))
            layer.CreateFeature(feature)
            #line.Destroy()
            #feature.Destroy()
        shapeData.Destroy()
    
    def __write_mesh_shp__(self, path):
        '''Wrtie the mesh to a shapefile'''
        #osgeo setup format and other stuff
        DriverFormat='ESRI Shapefile'    
        sr = osgeo.osr.SpatialReference()
        sr.ImportFromEPSG(int(self.wat.epsg))
        driver = osgeo.ogr.GetDriverByName(DriverFormat)
        if os.path.exists(path):
            driver.DeleteDataSource(path)
        shapeData = driver.CreateDataSource(path)
        layer = shapeData.CreateLayer('layer1', sr, osgeo.ogr.wkbPolygon)
        layerDefinition = layer.GetLayerDefn()
        #Fields and format to include in the shapefile
        names = ['polygon','area[m2]','z[m]','nfaces']
        fmts = [osgeo.ogr.OFTInteger,osgeo.ogr.OFTReal,osgeo.ogr.OFTReal,
              osgeo.ogr.OFTInteger]
        for name,fmt in zip(names, fmts):
            new_field=osgeo.ogr.FieldDefn(name,fmt)
            layer.CreateField(new_field)
        #Write polygon by polygon
        featureFID = 0
        for c, p in zip(self.polygons_id,self.polygons_topology):
            ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
            for i in np.array(p[-1]):
                ring.AddPoint(x=float(i[0]), y = float(i[1]))
            poly=osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feature = osgeo.ogr.Feature(layerDefinition)
            feature.SetGeometry(poly)
            feature.SetFID(0)
            feature.SetField('polygon',int(c))
            feature.SetField('area[m2]',float(p[2]))
            feature.SetField('z[m]',float(p[1]))
            feature.SetField('nfaces',int(p[3]))
            layer.CreateFeature(feature)
        shapeData.Destroy()
    
    def __get_left_right__(self, poly_xy, centers):
        left_right = []
        for c in centers:
            dist = np.linalg.norm(poly_xy - c, ord = 2, axis=1)  
            pos = []
            while len(pos)<2:
                p = np.argmin(dist)
                pos.append(p+1)
                dist[p] = 9999
            left_right.append(pos)
        self.river_left_right = left_right
    
    def __polygon_area__(self, x,y):
        correction = x[-1] * y[0] - y[-1]* x[0]
        main_area = np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:])
        return 0.5*np.abs(main_area + correction)

    def __get_z_from_dem__(self, x,y, kernel = 3, method = 'mean', percentile = 50):
        '''Get the elevation of a polygon based on the DEM'''
        #Get the col and the row that correspond to the polygon centriod
        col = np.floor((x - self.dem_prop[2])/self.dem_prop[-3]) -1 
        row = self.dem_prop[1] - np.floor((y-self.dem_prop[3])/self.dem_prop[-2])+1
        #Computes the elevation using different methods
        if method == 'mean':
            return np.mean(self.DEM[int(col)-kernel:int(col)+kernel, int(row)-kernel:int(row)+kernel])
        if method == 'percentile':
            return np.percentile(self.DEM[int(col)-kernel:int(col)+kernel, int(row)-kernel:int(row)+kernel], percentile)
        if method =='min':
            return np.min(self.DEM[int(col)-kernel:int(col)+kernel, int(row)-kernel:int(row)+kernel])
        if method =='max':
            return np.max(self.DEM[int(col)-kernel:int(col)+kernel, int(row)-kernel:int(row)+kernel])
    
    def __clean_river_points__(self, min_dist = 250):
        h_orders = np.array(self.river_topology)[:,-2]
        max_h = h_orders.max()
        drop_list = []
        excluded_pos = []
        for order in np.arange(max_h, 0, -1):
            pos = np.where(h_orders == order)[0]
            if order == max_h:
                pos = pos[1:]
            for c,i in enumerate(self.mesh_points_river.T[pos]):
                excluded_pos.append(pos[c])
                dist = np.linalg.norm(self.mesh_points_river.T - i, ord = 2, axis=1)    
                for c2,d in enumerate(dist):
                    if d>0 and d<min_dist and c2 not in excluded_pos:
                        drop_list.append(c2)
                        excluded_pos.append(c2)
        self.mesh_points_river = np.delete(self.mesh_points_river.T, drop_list, axis = 0).T
    
    def __clean_mesh_points__(self, min_dist = 100):        
        if self.focus_river is None:
            print('No focus definition, will proceed with %.2f distance for all the segments' % min_dist)
            pos = []        
            for i in self.mesh_points_river.T:
                dist = np.linalg.norm(self.mesh_points_dem.T - i, ord = 2, axis=1)    
                if np.min(dist) == 0:
                    dist[dist == 0] = 9999
                if np.min(dist) < min_dist:
                    pos.append(np.argmin(dist))
            self.mesh_points_dem = np.vstack([np.delete(self.mesh_points_dem[0], pos),np.delete(self.mesh_points_dem[1], pos)])
        else:
            print('Using focus setup, will use distance based on the config if avaiable')
            pos = []        
            for category, i in zip(self.focus_river, self.mesh_points_river.T):
                #Compute the distance 
                dist = np.linalg.norm(self.mesh_points_dem.T - i, ord = 2, axis=1)    
                if np.min(dist) == 0:
                    dist[dist == 0] = 9999
                try:
                    min_dist_temp= self.focus_dict[str(category)]['min_dem2river_distance']
                except:
                    min_dist_temp = min_dist
                if np.min(dist) < min_dist_temp:
                    pos.append(np.argmin(dist))
            self.mesh_points_dem = np.vstack([np.delete(self.mesh_points_dem[0], pos),np.delete(self.mesh_points_dem[1], pos)])
    
    def __get_segments_center_length__(self):
        '''Obtains the X,Y centroid and the straight length of each 
        segment.
        Parameters:
            - prop: list of list with the topology of the segments.
        Returns:
            - self.river_centers: list with the centers size: prop -1
            - self.river_length: list with the length size: prop -1'''
        centers = []
        lenght = []
        for l in range(1,len(self.river_topology)):
            xo = self.river_topology[l][2]
            yo = self.river_topology[l][3]
            d = self.river_topology[l][1]
            xd = self.river_topology[d][2]
            yd = self.river_topology[d][3]
            xc = (xo+xd)/2.
            yc = (yo+yd)/2.
            centers.append([xc,yc])
            lenght.append(np.sqrt((xo-xd)**2 + (yo-yd)**2))
        self.river_center = centers
        self.river_length = lenght
        
    def __get_segment_sinuosity__(self):
        '''Obtains the sinuosity of the segments using the straight length 
        end the cells directional path.
        Parameters: 
            - prop: list of list with the topology of the segments.
            - length: list with the segments length (derived from get_segments_length).
        Returns:
            - sinuosity: the sinuosity factor (Lr / Ls).
            - real_length: the lenght of the link crossing the DEM cells.'''
        sinuos = []
        real_length = []
        for p,le in zip(self.river_topology[1:], self.river_length):
            pos_o = p[-1]
            pos_d = self.river_topology[p[1]][-1]
            r_length = 0
            while pos_o != pos_d:
                r_length += self.wat.CellLong[pos_o]
                pos_o = self.wat.ncells - self.wat.structure[0][pos_o]
            sino = r_length / le
            if sino < 1:
                sino = 1
            sinuos.append(sino)
            real_length.append(r_length)
        self.river_sinuosity =  sinuos
        self.river_real_length = real_length
    
    def __channel2segments__(self, link, start = 1, dest = -4):
        '''Converts a channel link to a set of N segments 
        Parameters: 
            - link: the link to process.
            - start: the initial count for the segments.
            - dest: the downstream segment.            
        Returns:
            -segment list: id, dest, x_start, y_start, z_min, z_max, order, cell_start'''
        #Get the hydrological segment
        pos = np.where(self.links == link)
        if pos[0].size > 0:
            stream = self.wat.CellLong[pos]
        else:
            pos = np.where(self.links2 == link)
            stream = self.wat.CellLong[pos]

        #Determines the focus group at which the link belongs 
        if self.focus_map is not None:
            try:
                category = int(stats.mode(self.focus_map[pos]).mode[0])
                group = str(category)
                threshold = self.focus_dict[group]['seg_threshold']
            except:
                print('Warning: it is possible that %d is not defined in self.focus_dict or' % category)
                print(' that the property seg_threshold has not been defined')
                print('using the default self.threshold')
        else:
            threshold = self.threshold
            category = 0
        
        #Obtains the segments of the channel 
        stream_cat = np.ceil(stream.cumsum() /threshold)
        last = stream_cat[-1]
        stream_cat = last+1-stream_cat

        #Border correction of small channels 
        if stream_cat.size > 1:
            if stream_cat[stream_cat == stream_cat[-1]].size == 1:
                stream_cat[-1] = stream_cat[-2]
                stream_cat = stream_cat - 1 
            if stream_cat[stream_cat == stream_cat[0]].size == 1:
                stream_cat[0] = stream_cat[1]

        #Set the properties of the segment
        xt = self.x[pos]
        yt = self.y[pos]
        Z = self.wat.CellHeight[pos]
        H = self.wat.CellHorton_Stream[pos]

        properties = []
        focus_categories =[]
        for c, i in enumerate(np.unique(stream_cat)):
            pos2 = np.where(stream_cat == i)[0]
            prop = []
            prop = [int(i+start)]        
            if i > stream_cat[-1]:
                prop.append(int(i-1+start))
            else:
                prop.append(int(dest))
            #SEt the x, y coordinates of the link start
            prop.append(xt[pos2][0])
            prop.append(yt[pos2][0])
            #Get the height of the link 
            prop.append(Z[pos2].mean() - 20)
            prop.append(Z[pos2].mean())
            #Get the order 
            prop.append(int(np.max(H[pos2])))
            #Get the cell position of the start of the segment 
            prop.append(pos[0][pos2][0])
            properties.append(prop)
            #Update the focus category for that segment
            focus_categories.append(category)
            
        return properties, xt[pos2][0], yt[pos2][0], int(i+start), focus_categories
    
    def read_attfile(self, path):
        f = open(path,'r')
        lines = f.readlines()
        f.close()
        z = [i.split()[:8] for i in lines[1:]]
        df = pd.DataFrame(z, columns=lines[0].split(' ')[:8])
        df.set_index('INDEX', inplace = True)
        df.index = df.index.astype(int)
        df = df.astype(int)
        self.attr = df
        
    def __neighbors2array__(self):
        '''Obtains a list with the neigbors as arrays'''
        self.neighbors = []
        for p in self.polygons_topology:
            self.neighbors.append(np.array(p[4]) * np.array(p[7]))

        #Eliminates neighbors that are not in the self.polygons_topology list
        for c, n in enumerate(self.neighbors):
            if n.max() > self.polygons_id[-1]:
                for c2, i in enumerate(n):
                    if i > self.polygons_id[-1]:
                        self.neighbors[c][c2] = 0
                        
    def __correct_numbering__(self, verbose = False):
        '''Corrects the ids of the polygons taking away the ids that jump'''
        #Define the lists with the new polygons ids and their neighbors
        new_ids = []
        new_neighbors = copy.deepcopy(self.neighbors)
        c = 1
        c2 = 1
        #Progress bar
        print('Renumbering polygons ids ...')
        out = display(progress(0, len(self.polygons_id)), 
                    display_id=True)
        #Iterates through all the polygons
        for pid in self.polygons_id:
            #If the polygon id is greater than the expected ID, changes it in the neighbors
            if pid > c:
                #Get the new id 
                new_pid = c
                new_ids.append(new_pid)        
                #Prints if verbose
                if pid >c2 and verbose:
                    print(pid)
                    c2 = pid
                #Correct the number in the neighbors
                for n1, n2 in zip(self.neighbors, new_neighbors):
                    #Change the poly id if it is on that neighbor
                    if n1[n1==pid].size > 0:
                        n2[n2==pid] = new_pid
            else:
                new_ids.append(pid)
            out.update(progress(c, len(self.polygons_id)))
            c+=1
            c2+=1      
        #Replaces the polygons ids in self.polygon_ids and self.neighbors
        print('Done renumbering')
        self.neighbors = new_neighbors
        self.polygons_id = new_ids
                        
    def check_if_neighbors_exist(self):
        '''Checks if any polygon has a neighbor with an ID higher than the total number of elements'''
        #Set the number of elements and an empty list of bad polygons
        npolygons = len(self.polygons_topology)
        bad_neighbors = []
        #Iterate trhough all the polygons
        for c, polygon in enumerate(self.polygons_topology):
            neighbors = np.array(polygon[4])
            existing = np.array(polygon[7])
            neighbors = neighbors * existing
            if neighbors.max() > npolygons:
                bad_neighbors.append([c+1, neighbors.argmax()])
        #Returns a list that has the polygon followed by the polygon Id that is outside of the range
        return bad_neighbors
