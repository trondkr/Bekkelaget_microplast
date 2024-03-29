# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import logging
import numpy as np
from opendrift.models.oceandrift3D import OceanDrift3D
from opendrift.elements.passivetracer import PassiveTracer

# We add the property 'wind_drift_factor' to the element class
PassiveTracer.variables = PassiveTracer.add_variables([
    ('wind_drift_factor', {'dtype': np.float32,
                           'unit': '%',
                           'default': 0.02}),
    ('depth', {'dtype': np.float32,
               'unit': '%',
               'default': -16.0}),
    ('terminal_velocity', {'dtype': np.float32,
                           'units': 'm/s',
                           'default': 0.01})])


class PlastDrift(OceanDrift3D):
    """Trajectory model based on the OpenDrift framework.

    Propagation of plastics particles with ocean currents and
    additional Stokes drift and wind drag.

    Developed at MET Norway.

    """

    ElementType = PassiveTracer

    required_variables = [
        'x_sea_water_velocity', 'y_sea_water_velocity',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_significant_height',
        'sea_floor_depth_below_sea_level',
        'x_wind', 'y_wind',
        'ocean_vertical_diffusivity']
    required_variables.append('land_binary_mask')

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'sea_surface_wave_stokes_drift_x_velocity': 0,
                       'sea_surface_wave_stokes_drift_y_velocity': 0,
                       'x_wind': 0,
                       'y_wind': 0,
                       'sea_surface_wave_significant_height': 0,
                       'ocean_vertical_diffusivity': .02}

    required_profiles = None

    configspecPlastDrift = '''
        [general]
            coastline_action = option('none', 'stranding', 'previous', default='previous')
        [drift]
            verticaladvection = boolean(default=True)
            use_tabularised_stokes_drift = boolean(default=True)
            wind_drift_depth = float(min=0, max=10, default=0.1)
        [turbulentmixing]
            diffusivitymodel = option('environment', 'stepfunction', 'windspeed_Sundby1983', 'windspeed_Large1994', 'gls_tke', default='windspeed_Large1994')
        '''

    def __init__(self, *args, **kwargs):

        # Call parent constructor
        super(PlastDrift, self).__init__(*args, **kwargs)

        # Add specific config settings
        self._add_configstring(self.configspecPlastDrift)


    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        self.update_particle_depth_sinking()

        # Advect particles due to Stokes drift
        self.stokes_drift()

        # Advect particles due to wind-induced shear near surface
        self.advect_wind()

       
    def update_particle_depth_sinking(self):

        logging.debug('Submerging according to sinking velocity')
        self.elements.z = self.elements.z+self.elements.terminal_velocity*self.time_step.seconds
       # logging.debug("Particle 1 z: {}".format(self.elements.z[0]))
       # logging.debug("Z: {} vel: {} seconds: {}".format(self.elements.z,self.elements.terminal_velocity,self.time_step.seconds))

    def update_particle_depth(self):

        w = self.wind_speed()
        if w.max() == 0:
            logging.debug('No wind, all ellements at surface')
            self.elements.z = np.zeros(len(self.elements.z))
            return
        logging.debug('Submerging according to wind')
        self.elements.z = -np.random.exponential(
            scale=self.environment.ocean_vertical_diffusivity/
                    self.elements.terminal_velocity,
            size=self.num_elements_active())
