#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.models.oceandrift import OceanDrift
from IBM.pelagicplankton import PelagicPlanktonDrift
from pprint import pprint
from netCDF4 import Dataset, date2num, num2date
from scipy.ndimage.filters import gaussian_filter
import matplotlib
import os
import config_bekkelaget as confm
import probability_distribution_map_v2 as pb
import time

o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information

#######################
# Preparing readers
#######################

confobj=confm.bekkelaget_conf()

base = 'results'
baseout = 'figures'
hexagon = False
startdate = '01062015'
enddate = '30122015'
experiment = 1

start_time = time.time()
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
            confobj.plot_type='drift'

            reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=10.7, llcrnrlat=59.8,
                                urcrnrlon=10.8, urcrnrlat=59.9,resolution='h', 
                                projection='merc')

            o.add_reader(reader_basemap)

            plotfilenameAnime = 'figures/Bekkelaget_{}_drift_sinkspeed_{}_experiment_{}.mp4'.format(confobj.species,
                abs(confobj.select_sinking_velocity),
                confobj.experiment)
            plotfilenameColor = 'figures/Bekkelaget_{}_drift_sinkspeed_{}_experiment_{}_color.png'.format(confobj.species,
                abs(confobj.select_sinking_velocity),
                confobj.experiment)   
            plotfilename = 'figures/Bekkelaget_{}_drift_sinkspeed_{}_experiment_{}.png'.format(confobj.species,
                abs(confobj.select_sinking_velocity),
                confobj.experiment) 
            
            pb.createInputFilename(confobj)

            if os.path.exists(confobj.outputFilename):
                print("=> Opening input file: {}".format(confobj.outputFilename)
                )
                o.io_import_file(confobj.outputFilename)
                # o.add_reader([reader_basemap]) #Do not include basemap when stranding is deactivated

            # o.plot_vertical_distribution()
            #  o.plot(linecolor='z',lvmin=-150, lvmax=0,filename=plotfilenameColor)
                o.plot(filename=plotfilename)
                o.animation(filename=plotfilenameAnime)
            else:
                print("=> File does not exist {}".format(confobj.outputFilename))