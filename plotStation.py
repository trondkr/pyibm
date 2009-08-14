#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 8, 13)
__modified__ = datetime.datetime(2009, 8, 13)
__version__  = "0.1"
__status__   = "Development"

station="/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc_res.nc"
cdf=Dataset(station,"r")


light=cdf.variables["average_light"][:,:,:,:]
depth=cdf.variables["depth"][:,:,:,:]
time=(cdf.variables["time"][:])/24.
print time
sgr=cdf.variables["sgr"][:,:,:,:]
wgt=cdf.variables["wgt"][:,:,:,:]
length=cdf.variables["length"][:,:,:,:]
survival=cdf.variables["survival_probability"][:,:,:,:]

sgrrel=cdf.variables["sgr_rel"][:,:,:,:]

for ind in range(len(wgt[0,:,0,0])):
    plt.plot(time[0:400],wgt[0,ind,0:400,0], 'b-')



show()