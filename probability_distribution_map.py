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
import config_bekkelaget as confm
from fast_histogram import histogram2d
import cmocean
import laplacefilter
import mpl_util
from matplotlib.pyplot import cm 
import drawCircleOfDistance
import drawBathymetry

def createInputFilename(confobj):
    startDate=''
    if confobj.startdate.day<10:
        startDate+='0%s'%(confobj.startdate.day)
    else:
        startDate+='%s'%(confobj.startdate.day)

    if confobj.startdate.month<10:
        startDate+='0%s'%(confobj.startdate.month)
    else:
        startDate+='%s'%(confobj.startdate.month)

    startDate+='%s'%(confobj.startdate.year)

    endDate=''
    if confobj.enddate.day<10:
        endDate+='0%s'%(confobj.enddate.day)
    else:
        endDate+='%s'%(confobj.enddate.day)

    if confobj.enddate.month<10:
        endDate+='0%s'%(confobj.enddate.month)
    else:
        endDate+='%s'%(confobj.enddate.month)

    endDate+='%s'%(confobj.enddate.year)
 
    
    outputFilename='results/Bekkelaget_{}_drift_{}_to_{}_sinkspeed_{}_experiment_{}.nc'.format(confobj.species,
    startDate,endDate,
    abs(confobj.select_sinking_velocity),
    confobj.experiment)
    
    confobj.outputFilename=outputFilename

def get_paths(confobj):
    
    createInputFilename(confobj)
    return confobj.outputFilename

def create_map(confobj):
    plt.clf()
    confobj.ax = plt.subplot(111)
    mymap = Basemap(llcrnrlon=confobj.xmin, llcrnrlat=confobj.ymin,
                    urcrnrlon=confobj.xmax, urcrnrlat=confobj.ymax,resolution='h', 
                    projection='merc')

    mymap = Basemap(llcrnrlon=10.70, llcrnrlat=59.80,
                    urcrnrlon=10.80, urcrnrlat=59.9,resolution='f', 
                    projection='merc')

   # mymap.drawmapboundary(fill_color='#677a7a')
   
   # mymap.fillcontinents(color='lightgrey',zorder=2)
    mymap.drawcoastlines()
    confobj.mymap=mymap
  #  drawBathymetry.addBathymetry(confobj)

def get_pos(confobj):     
    print("Opening file {}".format(confobj.paths)) 
    df = xr.open_mfdataset(confobj.paths, concat_dim='trajectory') 
    d = df.groupby(df.trajectory) #.sel(time=slice(confobj.startdate, confobj.enddate))
  
    return df.lat[24:-1,:].values,df.lon[24:-1,:].values, df.z[24:-1,:].values, df.time[24:-1]

def distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295     #Pi/180
    a = 0.5 - np.cos((lat2 - lat1) * p)/2 + np.cos(lat1 * p) * np.cos(lat2 * p) * (1 - np.cos((lon2 - lon1) * p)) / 2
    return 12742 * np.arcsin(np.sqrt(a))

def createBins(confobj):

    print('func: createBins() => Creating bins for averaging')
   
    dy = distance(confobj.ymin,confobj.xmin,confobj.ymax,confobj.xmin)
    dx = distance(confobj.ymin,confobj.xmin,confobj.ymin,confobj.xmax)
    #dy = distance(59.8,10.7,59.9,10.7)
    #dx = distance(59.8,10.7,59.8,10.8)
    
    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))
    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))
    
    degrees_to_radians = np.pi/180.0
    radians_to_degrees = 180.0/np.pi
    earth_radius=6371
    # Distance in km between new bins converted to degrees
    
    r = earth_radius*np.cos(confobj.ymin*degrees_to_radians)
    confobj.deltaX = (confobj.requiredResolution/r)*radians_to_degrees
    confobj.lon_bins = np.arange(np.floor(confobj.xmin),np.ceil(confobj.xmax),0.01) 
    confobj.dx=len(confobj.lon_bins)
    print(confobj.lon_bins)

    confobj.deltaY = (confobj.requiredResolution/earth_radius)*radians_to_degrees
    confobj.lat_bins = np.arange(np.floor(confobj.ymin),np.ceil(confobj.ymax),0.01)
    #confobj.deltaY)
    confobj.dy=len(confobj.lat_bins)
    print(confobj.lat_bins)
    
    print('=> created binned array of domain of grid cell size (%s,%s) with resolution %s'%(confobj.deltaX,confobj.deltaY,confobj.requiredResolution))
    
    
def get_density(lats, lons, nlevels, confobj):

    XX=np.ma.masked_where(lons>200, lons)
    YY=np.ma.masked_where(lats>200, lats)
    print("test", XX.ravel().mask, len(XX.ravel()))
   # density,_,_ = np.histogram2d(XX.ravel(),YY.ravel(), bins=[confobj.dx,confobj.dy], density=False, weights=XX.ravel().mask)
    density = histogram2d(XX, YY, range=[[confobj.xmin,confobj.xmax],[confobj.ymin,confobj.ymax]], bins=[confobj.dx,confobj.dy])
    print(density)
    total = np.sum(density)
    print(len(lons),len(lats),np.shape(density))
    density=(density/total)*100.
    
    print("Total number of points {} percentage sum {}".format(total,np.sum(density)))
    print(np.min(density),np.max(density), np.mean(density))
   # density = ma.masked_where(density == 0, density)
    levels = MaxNLocator(nbins=nlevels).tick_values(0,2)
    levels=np.arange(density.min(),2,0.1)

    norm = BoundaryNorm(levels, ncolors=confobj.cmap.N, clip=True)

    # Turn the lon/lat of the bins into 2 dimensional arrays 
    lon_bins_2d, lat_bins_2d = np.meshgrid(confobj.lon_bins, confobj.lat_bins)
  
    return lon_bins_2d,lat_bins_2d,density,norm

def make_map(confobj):
    create_map(confobj)
    lats,lons,depths,times= get_pos(confobj)
    
    createBins(confobj)   
    confobj.results_startdate=times[0].values
    confobj.results_enddate=times[-1].values
    
    if confobj.plot_type == 'heatmap':
        nlevels = 100    
        confobj.cmap = plt.get_cmap(confobj.cmapname)     
     
        lon_bins_2d,lat_bins_2d,density,norm = get_density(lats[:,:], lons[:,:], nlevels, confobj)
        xs, ys = confobj.mymap(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh       
    
        cs = confobj.mymap.pcolormesh(lon_bins_2d, lat_bins_2d, density.T, latlon=True, cmap=confobj.cmap, norm=norm, edgecolors='face',alpha=0.8,linewidths=0.1)   

   
        # LOKI distribution/seed area
        X,Y = drawCircleOfDistance.createCircleAroundWithRadius(confobj.st_lats[0],confobj.st_lons[0],confobj.releaseRadius/1000.)
        confobj.mymap.plot(X,Y,latlon=True,marker=None,color='w',linewidth=0.2)

        # LOKI main point
        x,y = confobj.mymap(confobj.st_lons[0], confobj.st_lats[0])
        confobj.mymap.plot(x,y ,marker='D',color='w',markersize=0.1)
        plt.colorbar(cs, shrink=0.7)
        figname = r'MarMine_{}_{}_to_{}_heatmap.png'.format(confobj.species,confobj.results_startdate,confobj.results_enddate) 
     

    elif confobj.plot_type == 'scatter':
        x,y = confobj.mymap(lons,lats)
        confobj.mymap.scatter(x,y,alpha = 0.3,c = 'r',s = 1)

        # LOKI distribution/seed area
        X,Y = drawCircleOfDistance.createCircleAroundWithRadius(confobj.st_lats[0],confobj.st_lons[0],confobj.releaseRadius/1000.)
        confobj.mymap.plot(X,Y,latlon=True,marker=None,color='b',linewidth=0.2)

        # LOKI main point
        x,y = confobj.mymap(confobj.st_lons[0], confobj.st_lats[0])
        confobj.mymap.plot(x,y ,marker='D',color='b',markersize=0.2)
        plt.colorbar(cs, shrink=0.7)
     

        figname = r'MarMine_{}_{}_to_{}_scatter.png'.format(confobj.species,confobj.startdate,confobj.enddate) 
     
    plt.savefig(figname,format = 'png',dpi = 300)



def call_make_map(confobj):
    
    confobj.paths = get_paths(confobj) 
    make_map(confobj)

if __name__ == "__main__":
    start_time = time.time()

    confobj=confm.bekkelaget_conf()
    
    experiments = [1]
    years=[2018]

    for year in years:
        
        confobj.startdate=datetime(year,2,5,0,0,0)
        confobj.enddate=datetime(year,3,5,0,0,0)
        marine_organism="microplast"

        for confobj.select_sinking_velocity in confobj.sinkingvelocities:
            for experiment in experiments:

                confobj.experiment=experiment
                confobj.species=marine_organism
                confobj.plot_type='scatter'
                call_make_map(confobj) 
