from wmf import wmf 
import numpy as np 
from scipy.spatial import Voronoi, voronoi_plot_2d
from skimage.morphology import dilation, square
import os
import osgeo

class ghost_preprocess():
    
    def __init__(self, watershed, seg_threshold = 400, seg_point_distance = 100,thre_sens = 10):
        '''Defines the class to derive the topology required by ghost (.riv and .mesh files)
        Parameters:
            - watershed: A watershed element obtained with wmf.SimuBasin            
            - seg_threshold: minimum distance to divide channels into segments.
            - seg_point_distance: the distance of the perpendicular points to the channel segments
            - thre_sens: in meters, the sensitivity to transform from channels to segments
        Returns:
            - self.wat: a copy of the watershed element
            - self.x and self.y: the X and Y cooridnates of the watershed element
            - self.links and self.links2: element that identifies where channels are'''
        self.threshold = seg_threshold
        self.seg_point_distance = seg_point_distance
        self.wat = watershed
        #Get basic geomorphology 
        self.wat.GetGeo_Cell_Basics()
        #self.wat.GetGeo_StreamOrder()
        self.wat.GetGeo_StreamOrder(threshold=watershed.threshold)
        #Get the coordinates of each element and the links only 
        self.x, self.y = wmf.cu.basin_coordxy(self.wat.structure, self.wat.ncells)
        self.links = self.wat.CellCauce * self.wat.hills_own
        channels = np.zeros(cu.ncells)
        channels[self.wat.CellAcum>self.wat.threshold-thre_sens]=1
        self.links2 = channels * cu.hills_own
    
    def get_segments_topology(self):
        '''Using channel2segment obtains the topological connection 
        of the segments.
        Parameters: 
            - threshold: the minium req lenght (m) to obtain a segment.
        Returns:
            - self.river_topology: list of list describing the connections among segments.
            - self.river_centers: list with the centers size: prop -1
            - self.river_length: list with the length size: prop -1'''
        #Code starts
        start = 0
        prop = [[0,-999,self.x[-1],self.y[-1],-999,-999,self.wat.CellHorton_Hill[-1],self.wat.ncells]]
        new_dest = [0]
        for c, dest in enumerate(self.wat.hills[1][::-1]):
            d = new_dest[dest]
            p, x1, y2, start = self.__channel2segments__(c+1, start, d)
            prop.extend(p)
            new_dest.append(start)
        self.river_topology = prop
        self.__get_segments_center_length__()
        self.__get_segment_sinuosity__()
    
    def get_mesh_river_points(self, dist = 100):
        '''Obtains the two centers of the lines that are perpendicular 
        to each river segment'''
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
            #Compute the new X and Y
            x1 = xc - self.seg_point_distance/(np.sqrt(1+m**2))
            y1 = m*(x1-xc)+yc

            x2 = xc + self.seg_point_distance/(np.sqrt(1+m**2))
            y2 = m*(x2-xc)+yc

            Xp.append([x1,x2])
            Yp.append([y1,y2])
        Xv = np.array(Xp).T.reshape(len(Xp)*2)
        Yv = np.array(Yp).T.reshape(len(Yp)*2)
        XYr = np.vstack([Xv, Yv])
        self.mesh_points_river = XYr
    
    def get_mesh_grid_points(self, mesh_spaces = 8, border_iter = 3):
        '''Obtains the X and Y coordinates for the mesh and the borders of the mesh.
        Parameters:
            - mesh_spaces: step to take valuesin the grid given by the DEM.
            - border_iter: number of times to perform the dilation that generates the 
                border points
        Returns:
            - XYm: np.array (2,Nmesh) with the coordinates of the points inside the mesh.
            - XYb: np.array(2,Nborder) with the coordinates of the points in the border'''
        # Get the X and Y values fom the DEM mesh
        x_vect, y_vect = wmf.cu.basin_coordxy(self.wat.structure, self.wat.ncells)
        x_mask = wmf.cu.basin_float_var2map(self.wat.structure, x_vect, wmf.cu.ncols,wmf.cu.nrows, self.wat.ncells)
        y_mask = wmf.cu.basin_float_var2map(self.wat.structure, y_vect, wmf.cu.ncols,wmf.cu.nrows, self.wat.ncells)

        #Obtain a set of points spaced in the grid
        x_steps = np.arange(0,x_mask.shape[0],mesh_spaces)
        y_steps = np.arange(0,x_mask.shape[1],mesh_spaces)
        xv,yv = np.meshgrid(x_steps, y_steps)

        # Get the xy points inside the 
        XYm = []
        for i in x_steps:
            for j in y_steps:
                if x_mask[i,j] > 0:
                    XYm.append([x_mask[i,j],y_mask[i,j]])
        XYm = np.array(XYm).T

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
        return XYm, borders
    
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

        #Obtains the segments of the channel 
        stream_cat = np.ceil(stream.cumsum() / self.threshold)
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
        return properties, xt[pos2][0], yt[pos2][0], int(i+start)