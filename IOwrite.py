
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

def writeStationFile(grdSTATION,ntime,d,p,outfilename):
        
     if grdSTATION.Initialized is False:
        
        grdSTATION.Initialized=True
        
        if os.path.exists(outfilename):
            os.remove(outfilename)
        
        #f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
        f1 = Dataset(outfilename, mode='w', format='NETCDF4')
        f1.description="This is a IBM result file for a given station"
        f1.history = 'Created ' + time.ctime(time.time())
        f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
        f1.type='NetCDF4 classic created using pyibm' 
           
        # Define dimensions
        f1.createDimension('depth', grdSTATION.Nlarva)
        f1.createDimension('time', None)
        f1.createDimension('prey', grdSTATION.Nprey)

        d=[i for i in range(grdSTATION.Nlarva)]
        
        vnc=f1.createVariable('depth','d',('depth',),zlib=True)
        vnc.long_name = "depth" ;
        vnc.units = "meter"
        vnc[:]=d
       
        v_time = f1.createVariable('time', 'd', ('time',),zlib=True)
        v_time.long_name = 'Days since 01-01 00:00:00'
        v_time.units = 'days'
        v_time.field = 'time, scalar, series'
        v_time.calendar='standard'
        
        v=f1.createVariable('wgt', 'f', ('time','depth','prey'),zlib=True)
        v.long_name = "Weight at day of year , depth, given prey"
        v.units = "mg"
        
        v=f1.createVariable('length', 'f', ('time','depth','prey'),zlib=True)
        v.long_name = "Length at day of year , depth, given prey"
        v.units = "mm"
        
        v=f1.createVariable('sgr', 'f', ('time','depth','prey'),zlib=True)
        v.long_name = "SGR at day of year , depth, given prey"
        v.units = "%day-1"
        
        v=f1.createVariable('sgrAF', 'f', ('time','depth'),zlib=True)
        v.long_name = "Maximum growth rate at day of year , depth, given prey"
        v.units = "%day-1"
        
        v=f1.createVariable('sgr relative', 'f', ('time','depth','prey'),zlib=True)
        v.long_name = "Growth rate relative to maximum growth rate at day of year , depth, given prey"
        v.units = "%day-1"
        
        v=f1.createVariable('survival probability', 'f', ('time','depth','prey'),zlib=True)
        v.long_name = "Survival probablity"
        v.units = "day-1"
        
        v=f1.createVariable('average light', 'f', ('time','depth'),zlib=True)
        v.long_name = "average daily light value at day of year , depth"
        v.units = ""
        
        v=f1.createVariable('temp', 'f', ('time','depth'),zlib=True)
        v.long_name = "Temperature"
        v.units = "degrees Celsius"
        
     else:
        f1 = Dataset(outfilename, mode='a', format='NETCDF4', zlib=True)
    
     f1.variables['time'][ntime]   = grdSTATION.larvaTime 
     
     f1.variables['wgt'][ntime,d,p]  = grdSTATION.larvaWgt
     f1.variables['length'][ntime,d,p]  = grdSTATION.larvaLength
     f1.variables['sgr'][ntime,d,p]  = grdSTATION.larvaSgr
     f1.variables['average light'][ntime,d]  = grdSTATION.larvaAveLight
     f1.variables['sgrAF'][ntime,d]  = grdSTATION.larvaSgrAF
     f1.variables['temp'][ntime,d]  = grdSTATION.larvaTdata
     f1.variables['sgr relative'][ntime,d,p]  = (grdSTATION.larvaSgr/grdSTATION.larvaSgrAF)*100.
     f1.variables['survival probability'][ntime,d,p]  = grdSTATION.larvaPsur
     f1.close()


