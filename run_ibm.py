from datetime import datetime, timedelta
import numpy as np
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.plastdrift_bekkelaget import PlastDrift
import logging
import os
from netCDF4 import Dataset, date2num,num2date
from numpy.random import RandomState
import matplotlib.pyplot as plt
import time
import config_bekkelaget as confm
import random
import logging

__author__ = 'Trond Kristiansen'
__email__ = 'me (at) trondkristiansen.com'
__created__ = datetime(2019, 7, 26)
__modified__ = datetime(2019, 7, 26)
__version__ = "1.0"
__status__ = "Development, modified on 25.06.2018, 06.08.2018, 26.07.2019"

def createOutputFilenames(confobj):
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
    
    if not os.path.exists('results'):
        os.makedirs('results')
    if os.path.exists(outputFilename):
        os.remove(outputFilename)
        
    confobj.outputFilename=outputFilename

   
def createAndRunSimulation(confobj):

    # Setup a new simulation
    o = PlastDrift(loglevel=0)  # Set loglevel to 0 for debug information
    o.complexIBM=False

  #  reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=confobj.xmin, llcrnrlat=confobj.ymin,
  #                  urcrnrlon=confobj.xmax, urcrnrlat=confobj.ymax,resolution='f', 
  #                  projection='merc')
                       
 #  o.add_reader([reader_basemap]) 
   
    reader_physics = reader_ROMS_native.Reader(confobj.basedir+confobj.pattern)
    o.add_reader([reader_physics ]) 

    #######################
    #Adjusting configuration
    #######################
    o.set_config('processes:turbulentmixing',  True)
    o.set_config('processes:verticaladvection', True)
    o.set_config('turbulentmixing:diffusivitymodel','gls_tke')
    o.set_config('turbulentmixing:TSprofiles', True)
    #o.set_config('turbulentmixing:mixingmodel','randomwalk')
    o.set_config('drift:scheme', 'runge-kutta4')
    o.set_config('general:coastline_action', 'previous') #Prevent stranding, jump back to previous position
    o.set_config('general:basemap_resolution', 'f')
    o.set_config('drift:lift_to_seafloor', True)
   
    stdev=(confobj.low_depth-confobj.mean_depth)/3.
    z_levels=[]
    while len(z_levels) < confobj.releaseParticles:
        sample=random.gauss(confobj.mean_depth, stdev)
        if sample >= confobj.low_depth and sample < confobj.high_depth:
            z_levels.append(sample)
        else:
            z_levels.append(confobj.mean_depth)

    print('Seeding {} elements within a radius of {} m between depths {} to {} m'.format(confobj.releaseParticles, 
    confobj.releaseRadius,
    np.min(z_levels),
    np.max(z_levels)))
   
    print("Releasing {} larvae between {} and {}".format(confobj.species,confobj.startdate,confobj.enddate))
    o.seed_elements(lon=confobj.st_lons, 
                    lat=confobj.st_lats, 
                    number=confobj.releaseParticles, 
                    radius=[confobj.releaseRadius], cone=False, 
                    time=[confobj.startdate,confobj.enddate],
                    terminal_velocity=confobj.sinkingvelocities[confobj.experiment],
                    z=z_levels)
    
    print('Elements scheduled for {} : {}'.format(confobj.species,o.elements_scheduled))

    o.run(end_time=confobj.enddate, 
          time_step=timedelta(minutes=30),
          time_step_output=timedelta(minutes=30),
          outfile=confobj.outputFilename,
          export_variables=['sea_floor_depth_below_sea_level','z','terminal_velocity'])


if __name__ == "__main__":
    start_time = time.time()
    
    experiments = [0,1] # Used as index other places so have to run from 0---N
    years=[2019]

    for experiment in experiments:

        confobj=confm.bekkelaget_conf(experiment)

        for year in years:
            
            confobj.startdate=datetime(year,7,1,12,0,0)
            confobj.enddate=datetime(year,7,7,12,0,0)
            logging.debug("Running experiment for period {} to {}".format(confobj.startdate.year,confobj.enddate))

            for confobj.select_sinking_velocity in confobj.sinkingvelocities:
                confobj.experiment=experiment
                createOutputFilenames(confobj)
                logging.debug("Result files will be stored as:\nnetCDF=> {}".format(confobj.outputFilename))
                createAndRunSimulation(confobj)

    logging.debug("---  It took %s seconds to run the script ---" % (time.time() - start_time)) 