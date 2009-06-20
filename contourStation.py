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

station="41.5423_-132.9234_res.nc"

cdf=Dataset(station,"r")

#lon=cdf.variables["lon"][:]
light=flipud(rot90(cdf.variables["average light"][:,:]))
depth=cdf.variables["depth"][:]
time=cdf.variables["time"][:]
sgr=flipud(rot90(cdf.variables["sgr"][:,:,0]))
wgt=flipud(rot90(cdf.variables["wgt"][:,:,0]))
survival=flipud(rot90(cdf.variables["survival probability"][:,:,0]))
print survival
sgrrel=flipud(rot90(cdf.variables["sgr relative"][:,:,0]))
Z=survival
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

CS2 = contour(X, Y, Z, CS.levels[::1],
                        colors = 'k',
                        origin=origin,
                        hold='on')

#title('Climatologial temperature at station (%2.1f, %2.1f))'%(lon,lat))
xlabel('Day of year')
ylabel('Depth (m)')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = colorbar(CS)
cbar.ax.set_ylabel('verbosity coefficient')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)


show()