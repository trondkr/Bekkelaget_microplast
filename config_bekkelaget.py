#
#
# Config object for MarMine
import time, calendar
from datetime import datetime, timedelta

__author__ = 'Trond Kristiansen'
__email__ = 'me (at) trondkristiansen.com'
__created__ = datetime(2019, 7, 31)
__modified__ = datetime(2019, 7, 31)
__version__ = "1.0"
__status__ = "Development, modified on 31.07.2019"


class bekkelaget_conf():

    def user_defined_inits(self,experiment):
       
        self.requiredResolution=0.1 # km between bins
        self.xmin=10.475
        self.xmax=10.85
        self.ymin=59.6
        self.ymax=59.95

        self.probxmin=10.70
        self.probxmax=10.78
        self.probymin=59.87
        self.probymax=59.90

        self.releaseParticles=500
        self.releaseRadius=50

        # Spring depth where outlet water is in balance with fjord water - mean 18 m
        if experiment==0:
            self.low_depth=-35
            self.high_depth=-1
            self.mean_depth= -18.0

        # Fall - mean 32 m
        if experiment==1:
            self.low_depth=-47.5
            self.high_depth=-16.5
            self.mean_depth= -32.0

        self.verticalBehavior=False
        self.basedir='/cluster/projects/nn9297k/FjordOS/' 
        self.pattern='ocean_his.nc_*' 
        self.species='microplast'
        self.plot_type='scatter'
        self.cmapname='RdYlBu_r'
        self.selectyear='all'
        # Bekkelaget - seed locations
        self.st_lons = [10.7528]
        self.st_lats = [59.8830]
        self.sinkingvelocities=[-0.0005, -0.005, -0.01] #m/s for three sizes of microplastic
        self.plastlengths=[0.0005,0.001,0.002] # mm
        self.interpolation='linearNDFast'
        self.sed_crit = 0.1

    def define_orgamisms(self):
        return {'plast': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']}

    def __init__(self,experiment):
        print('\n--------------------------\n')
        print('Started ' + time.ctime(time.time()))

        self.user_defined_inits(experiment)
        
        self.paths=None
        self.mymap=None
        self.ax=None
        self.deltaX = None
        self.deltaY = None
        self.dx=None
        self.dy=None
        self.cmap=None
        self.outputFilename=None
        self.results_startdate=None
        self.results_enddate=None
        self.startdate=None
        self.enddate=None   
        self.select_sinking_velocity=None