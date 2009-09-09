import netCDF4
from netCDF4 import Dataset
import types, math
import numpy as np
import datetime as dt
__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = dt.datetime(2008, 6, 10)
__modified__ = dt.datetime(2009, 6, 19)
__version__  = "1.2"
__status__   = "Production"

def getStationData(cdf, varlist, grdSTATION, log, clim):
    """
    This routine reads a netCDF4 file created using IOstation.py in soda2roms. Such
    a file only contain information for one lat/long location and is therefore different
    from finding values in a grid. Trond Kristiansen, 4.6.2009 on flight to Raleigh/Durham CO330
    """
    Nvars=len(varlist)
    t= (grdSTATION.endIndex-grdSTATION.startIndex)
    var_array_rawXY=np.zeros((t, Nvars),dtype=np.float64)
    var_array_rawXYZ=np.zeros((t, int(len(grdSTATION.depth)), Nvars),dtype=np.float64)
    
    var_numberXYZ=0
    var_numberXY =0
    
    if clim is False:
        grdSTATION.time=cdf.variables['time'][grdSTATION.startIndex:grdSTATION.endIndex]
    else:
        grdSTATION.time=(cdf.variables['clim_time'][grdSTATION.startIndex:grdSTATION.endIndex])*5
        
    for var in varlist:
        if var in ["taux","tauy","ssh"]:
            var_array_rawXY[:,var_numberXY] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex]
            var_numberXY=var_numberXY+1
        else:
            var_array_rawXYZ[:,:,var_numberXYZ] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex,:]
            var_numberXYZ=var_numberXYZ+1
   
        """
        Slice the data in the netcdf file that you want to keep. In our case we store a time-series from
        a station (eta_rho, xi_rho), from top to bottom. The variables extracted at this station are
        defined in the "varlist" variable.
        """
        
        if log is True:
            print "\n---> Extracting time series of %s from station"%(var)
            if var in ["taux","tauy","ssh"]:
                print '---> Maximum %s %f'%(var, var_array_rawXY[:,int(var_numberXY-1)].max())
                print '---> Minimum %s %f'%(var, var_array_rawXY[:,int(var_numberXY-1)].min())    
            else:
                print '---> Maximum %s %f'%(var, var_array_rawXYZ[:,:,int(var_numberXYZ-1)].max())
                print '---> Minimum %s %f'%(var, var_array_rawXYZ[:,:,int(var_numberXYZ-1)].min())      
   
    grdSTATION.data=var_array_rawXYZ
    grdSTATION.dataXY=var_array_rawXY
   