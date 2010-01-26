#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 6, 18)
__modified__ = datetime.datetime(2009, 6, 18)
__version__  = "0.1"
__status__   = "Development"

station="41.6423_-67.2001_res.nc"
station="/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc_res.nc"
#station="43.4111_-50.4321_res.nc"
cdf=Dataset(station,"r")

#lon=cdf.variables["lon"][:]
light=flipud(rot90(cdf.variables["average_light"][:,:,:,:]))
depth=cdf.variables["depth"][:,:,:,:]
time=cdf.variables["time"][:]
sgr=flipud(rot90(cdf.variables["sgr"][:,:,:,:]))
wgt=flipud(rot90(cdf.variables["wgt"][:,:,:,:]))
length=flipud(rot90(cdf.variables["length"][:,:,:,:]))
survival=flipud(rot90(cdf.variables["survival_probability"][:,:,:,:]))

sgrrel=flipud(rot90(cdf.variables["sgr_rel"][:,:,:,:]))
Z=length[:,:,0]
var='length'

origin = 'lower'

X=time
Y=-depth
    
CS = contourf(X, Y, Z,  10,
                        alpha=1.0,
                        cmap=cm.jet,
                        origin=origin)

# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution.

CS2 = contour(X, Y, Z, CS.levels,
                        colors = 'k',
                        origin=origin,
                        hold='on')

#title('Climatologial temperature at station (%2.1f, %2.1f))'%(lon,lat))
xlabel('Day of year')
ylabel('Depth (m)')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = colorbar(CS)
cbar.ax.set_ylabel('%s'%(var))
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)
plt.axis("tight")

show()