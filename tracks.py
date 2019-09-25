import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.pyplot import cm 
from mpl_toolkits.basemap import Basemap
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import animateScatter
import laplacefilter
import mpl_util
class Tracks(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self,lons,lats,z,figname):
     
        self.z = z
        self.figname=figname
        self.llcrnrlon=3.0
        self.llcrnrlat=69.4
        self.urcrnrlon=23.0
        self.urcrnrlat=75

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.mymap = Basemap(llcrnrlon=self.llcrnrlon, llcrnrlat=self.llcrnrlat,
                    urcrnrlon=self.urcrnrlon, urcrnrlat=self.urcrnrlat,resolution='h', 
                    projection='merc')
        self.mymap.shadedrelief()
        self.mymap.fillcontinents(color='#b8a67d',zorder=2)
        self.mymap.drawcoastlines()

        self.addBathymetry()
        animateScatter.KelpPolygons(self.mymap,self.ax)
       
        self.x,self.y=self.mymap(lons,lats)

    def plot_tracks(self):

        for traj_index in range(len(self.x[:,0])):
            print("Plotting trajectory {}".format(traj_index))
            self.mymap.plot(self.x[traj_index,:],self.y[traj_index],c='k',alpha=0.5,linewidth=0.4)
        
        plt.savefig(self.figname,format = 'png',dpi = 300)

    def addBathymetry(self):

        etopo1='/Users/trondkr/Dropbox/NIVA/ETOPO1/ETOPO1_Ice_g_gmt4.grd'
        etopo=xr.open_dataset(etopo1)
        lats=etopo['y'].values
        lons=etopo['x'].values

        res = self.findSubsetIndices(self.llcrnrlat-5,self.urcrnrlat+5,self.llcrnrlon-40,self.urcrnrlon+10,lats,lons)
        res = [ int(x) for x in res ]
     
        bathy = etopo["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])].values
        bathySmoothed = laplacefilter.laplace_filter(bathy,M=None)
        lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]]) 
        
        levels=[-6000,-5000,-4000,-3500,-3000,-2500,-2000,-1750,-1500,-1250,-1000,-750,-500, -400, -300, -250, -200, -150, -100, -50, -10]
        x, y = self.mymap(lon,lat)
        etopo_contour = self.mymap.contourf(x,y,bathySmoothed,levels,
                       cmap=mpl_util.LevelColormap(levels,cmap=cm.Blues_r),
                       extend='upper',
                       alpha=1.0,
                       origin='lower')
    
    
    def findSubsetIndices(self,min_lat,max_lat,min_lon,max_lon,lats,lons):

        """Array to store the results returned from the function"""
        res=np.zeros((4),dtype=np.float64)
        minLon=min_lon; maxLon=max_lon

        distances1 = []; distances2 = []
        indices=[]; index=1

        for point in lats:
            s1 = max_lat-point # (vector subtract)
            s2 = min_lat-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        distances1 = []; distances2 = []; index=1

        for point in lons:
            s1 = maxLon-point # (vector subtract)
            s2 = minLon-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
        minJ=indices[1][2]
        maxJ=indices[0][2]
        minI=indices[3][2]
        maxI=indices[2][2]

        res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
        return res
