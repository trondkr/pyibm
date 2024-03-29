
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

def writeStationFile(deltaH,deltaZ,grdSTATION,larvaTime,W,L,SGR,larvaTdata,larvaNauplii,larvaDepth,W_AF,larvaPsur,outputFile,startAndStop,Tlarva,Mfish,Minve,Mstar):

     f1 = Dataset(outputFile, mode='w', format='NETCDF4')
     f1.description="This is a IBM result file for a given station"
     f1.history = 'Created ' + time.ctime(time.time())
     f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
     f1.type='NetCDF4 classic created using pyibm'

     myoptions=""
     for opts in grdSTATION.OPTIONS:
        myoptions=myoptions + "%s "%(str(opts).rstrip())
     f1.options='IBM run with options: %s'%(myoptions)

     # Define dimensions
     f1.createDimension('cohort',grdSTATION.Ncohorts)
     f1.createDimension('larva', grdSTATION.Nlarva)
     f1.createDimension('time', None)
     f1.createDimension('larvaTime', Tlarva)
     f1.createDimension('IO', 2)

     v_time = f1.createVariable('time', 'd', ('time',),zlib=True)
     v_time.long_name = 'Seconds since 1948-01-01 00:00:00'
     v_time.units = 'seconds'
     v_time.field = 'time, scalar, series'
     v_time.calendar='standard'

     v=f1.createVariable('deltaH','d')
     v.long_name = "Hours per timestep" ;
     v.units = "1.0/hour"
     v[:]=deltaH

     v=f1.createVariable('deltaZ','d')
     v.long_name = "Vertical increments allowed for behavior" ;
     v.units = "meter"
     v[:]=deltaZ

     v=f1.createVariable('wgt', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Weight at day of year and depth"
     v.units = "mg"

     v=f1.createVariable('length', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Length at day of year , depth, given prey"
     v.units = "mm"

     v=f1.createVariable('sgr', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "SGR at day of year and depth"
     v.units = "% dt-1"

     v=f1.createVariable('sgrAF', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Maximum growth rate at day of year and depth"
     v.units = "%dt-1"

     v=f1.createVariable('sgr_rel', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Growth rate relative to maximum growth rate at day of year and depth"
     v.units = "%dt-1"

     v=f1.createVariable('survival_probability', 'f',('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Survival probablity"
     v.units = "dt-1"

     v=f1.createVariable('Mfish', 'f',('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Mortality from piscivores"
     v.units = "dt-1"

     v=f1.createVariable('Minvertebrates', 'f',('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Mortality from invertebrates"
     v.units = "dt-1"

     v=f1.createVariable('Mstarvation', 'f',('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Mortality from starvation"
     v.units = "dt-1"

     v=f1.createVariable('temp', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Temperature"
     v.units = "degrees Celsius"

     v=f1.createVariable('nauplii', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Nauplii"
     v.units = "numbers per liter"

     v=f1.createVariable('depth', 'f', ('cohort','larva','larvaTime'),zlib=True)
     v.long_name = "Depth"
     v.units = "meter"

     v=f1.createVariable('timeIndex', 'f', ('cohort','larva','IO'),zlib=True)
     v.long_name = "Array that gives the time index for entrance and exit in time of larvae"
     v.units = "index"

     f1.variables['time'][:]   = larvaTime

     f1.variables['wgt'][:,:,:]  = W[:,:,:]
     f1.variables['length'][:,:,:]  = L[:,:,:]
     f1.variables['sgr'][:,:,:]  = SGR[:,:,:]
     f1.variables['sgrAF'][:,:,:]  = W_AF
     f1.variables['temp'][:,:,:]  = larvaTdata
     f1.variables['nauplii'][:,:,:]  = larvaNauplii[:,:,:]
     f1.variables['depth'][:,:,:]  = larvaDepth
     f1.variables['timeIndex'][:,:,:]  = startAndStop

    # f1.variables['sgr_rel'][:,:,:,:]  = SGR/grdSTATION.larvaSgrAF)*100.
     f1.variables['survival_probability'][:,:,:]  = larvaPsur
     f1.variables['Mfish'][:,:,:]  = Mfish
     f1.variables['Minvertebrates'][:,:,:]  = Minve
     f1.variables['Mstarvation'][:,:,:]  = Mstar
     f1.close()
