# coding=utf-8

import os, sys
import numpy as np
import numpy.ma as ma
import glob
import matplotlib
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import xarray as xr
from datetime import datetime
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import animateScatter
import time
import utils


def find_depth(data):
    data['z'] = data['z'] * -1.
    # ! find first non nan at first and cut the rest 
    #data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    # find differences between floor depth and particle depth for each trajectory
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data

def get_groups(new_df,p_part):
    d = new_df.where(new_df.plantpart == p_part,drop = True)
    # apply method to each trajectory (particle release event)
    return d.groupby(d.trajectory).apply(find_depth)

def create_map():
    fig = plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,resolution='h', 
                    projection='merc')

    #mymap.drawmapboundary(fill_color='#677a7a')
    mymap.fillcontinents(color='#b8a67d',zorder=2)
    mymap.drawcoastlines()
    return mymap, ax

def get_pos(paths,kelpType):      
    df = xr.open_mfdataset(paths,concat_dim='trajectory')
    if kelpType != 'All':
        df = df.where(df.plantpart == kelpType,drop = True)     
    df = df.where(df.status > -1, drop = True)
    d = df.groupby(df.trajectory).apply(find_depth)
    parts = range(0,len(d.trajectory)-1)
    lats = [utils.get_lat(d,n).values for n in parts if utils.is_sedimented(d,n)]
    lons = [utils.get_lon(d,n).values for n in parts if utils.is_sedimented(d,n)]        
    return lats,lons

def distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295     #Pi/180
    a = 0.5 - np.cos((lat2 - lat1) * p)/2 + np.cos(lat1 * p) * np.cos(lat2 * p) * (1 - np.cos((lon2 - lon1) * p)) / 2
    return 12742 * np.arcsin(np.sqrt(a))

def createBins(requiredResolution=0.25):

    print('func: createBins() => Creating bins for averaging')
    xmin=-25.0; xmax=21.0
    ymin=60.0; ymax=76.0
   
    dy = distance(ymin,xmin,ymax,xmin)
    dx = distance(ymin,xmin,ymin,xmax)

    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))
    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))
    
    delta = int(np.round(dx/requiredResolution,0))
    print("delta {}".format(delta))
    xi = np.linspace(np.floor(xmin),np.ceil(xmax),delta)
    deltaX=abs(xi[0]-xi[1])
    print(deltaX)

    # We are only using longitude to calculate the approximate distance in degrees 
    # that is equvalent to requiredResolsution. Then we create latitude grid based on the 
    # same deltaX.
    yi = np.arange(np.floor(ymin),np.ceil(ymax),deltaX)

    print('=> created binned array of domain of grid cell size (%s) with resolution %s'%(delta,requiredResolution))
    
    return xi,yi,deltaX

#def get_bins(l,nbins):
#    db = 1.e-6 # bin padding    
#    return np.linspace(min(l)-db, max(l)+db, nbins)


def get_density(lats, lons,nlevels,cmap):

    requiredResolution=1
    lon_bins,lat_bins,deltaDegrees = createBins(requiredResolution)   

    density, _, _ = np.histogram2d(lats, lons, [lat_bins, lon_bins])
    density = ma.masked_where(density == 0, density)

    levels = MaxNLocator(nbins=nlevels).tick_values(0,100)
    #norm .tick_values(density.min(),density.max())
   
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # Turn the lon/lat of the bins into 2 dimensional arrays 
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
  
    return lon_bins_2d,lat_bins_2d,density,norm

def make_map(paths,kelpType,type,experiment,polygons):
    mymap, ax = create_map()
    lats,lons = get_pos(paths,kelpType)

    if type == 'heatmap':
        nlevels = 50        
        cmap = plt.get_cmap('Spectral_r')        
        lon_bins_2d,lat_bins_2d,density,norm = get_density(lats, lons,nlevels,cmap)
        xs, ys = mymap(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh          
        cs = mymap.pcolormesh(xs, ys, density, cmap=cmap, norm=norm)     
        polygons = animateScatter.KelpPolygons(mymap,ax,fill=False)
 
        plt.colorbar(cs, shrink=0.7)
        figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,str(experiment)) 
     

    elif type == 'scatter':
        x,y = mymap(lons,lats)
        mymap.scatter(x,y,alpha = 0.5,c = 'k',s = 10)
        figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,experiment)         

    plt.savefig(figname,format = 'png',dpi = 300)



def call_make_map(kelpType,plot_type,experiments,polygons):
    if len(experiments) > 1:
        p = [utils.get_paths(polygons,exp)for exp in experiments]
        paths = [val for sublist in p for val in sublist]
    else:
        paths = utils.get_paths(experiment = experiments[0],polygons = polygons) 

    make_map(paths,kelpType,plot_type,experiment = experiments,polygons = polygons)

if __name__ == "__main__":
    start_time = time.time()

    #Examples: 
    # 1 ) 
    # takes  20 seconds   
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = [1], polygons = [1]) # 'All')
    alltypes = (0,1,2,4)  
    experiments = (1,2,3,4,5)
    # 2)
    #Will create map for Kelp 1, polygon 1, all experiment 
    for experiment in experiments:
        for kelpType in alltypes:
            call_make_map(kelpType = kelpType,plot_type = 'heatmap',experiments = [experiment], polygons = 'All') # 'All')

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment
    # takes 20 minutes to plot     
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   
    #call_make_map(kelpType = 2,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment    
    #call_make_map(kelpType = 'All',plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')
   
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   
