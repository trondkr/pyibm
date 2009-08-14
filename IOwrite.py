
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

def writeStationFile(grdSTATION,larvaTime,W,L,SGR,larvaTdata,W_AF,aveLight,larvaPsur,outputFile):
     print outputFile
     f1 = Dataset(outputFile, mode='w', format='NETCDF4')
     f1.description="This is a IBM result file for a given station"
     f1.history = 'Created ' + time.ctime(time.time())
     f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
     f1.type='NetCDF4 classic created using pyibm' 
        
     # Define dimensions
     f1.createDimension('cohort',grdSTATION.Ncohorts)
     f1.createDimension('larva', grdSTATION.Nlarva)
     f1.createDimension('time', None)
     f1.createDimension('prey', grdSTATION.Nprey)
    
     v_time = f1.createVariable('time', 'd', ('time',),zlib=True)
     v_time.long_name = 'Days since 1948-01-01 00:00:00'
     v_time.units = 'days'
     v_time.field = 'time, scalar, series'
     v_time.calendar='standard'
     
     v=f1.createVariable('wgt', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Weight at day of year , depth, given prey"
     v.units = "mg"
     
     v=f1.createVariable('length', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Length at day of year , depth, given prey"
     v.units = "mm"
     
     v=f1.createVariable('sgr', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "SGR at day of year , depth, given prey"
     v.units = "%day-1"
     
     v=f1.createVariable('sgrAF', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Maximum growth rate at day of year , depth, given prey"
     v.units = "%day-1"
     
     v=f1.createVariable('sgr_rel', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Growth rate relative to maximum growth rate at day of year , depth, given prey"
     v.units = "%day-1"
     
     v=f1.createVariable('survival_probability', 'f',('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Survival probablity"
     v.units = "day-1"
     
     v=f1.createVariable('average_light', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "average daily light value at day of year , depth"
     v.units = ""
     
     v=f1.createVariable('temp', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Temperature"
     v.units = "degrees Celsius"
     
     v=f1.createVariable('depth', 'f', ('cohort','larva','time','prey'),zlib=True)
     v.long_name = "Depth"
     v.units = "meter"
    
     f1.variables['time'][:]   = larvaTime 
     
     f1.variables['wgt'][:,:,:,:]  = W[:,:,:,:]
     f1.variables['length'][:,:,:,:]  = L[:,:,:,:]
     f1.variables['sgr'][:,:,:,:]  = SGR[:,:,:,:]
     f1.variables['average_light'][:,:,:,:]  = aveLight
     f1.variables['sgrAF'][:,:,:,:]  = W_AF
     f1.variables['temp'][:,:,:,:]  = larvaTdata
    # f1.variables['sgr_rel'][:,:,:,:]  = SGR/grdSTATION.larvaSgrAF)*100.
     f1.variables['survival_probability'][:,:,:,:]  = larvaPsur
     f1.close()


