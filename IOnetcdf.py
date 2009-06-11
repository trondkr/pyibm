import netCDF4
from netCDF4 import Dataset
import datetime, types, math
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 10)
__modified__ = datetime.datetime(2009, 6, 4)
__version__  = "1.2"
__status__   = "Production"

def getStationData(cdf, varlist, grdSTATION, log):
    """
    This routine reads a netCDF4 file created using IOstation.py in soda2roms. Such
    a file only contain information for one lat/long location and is therefore different
    from finding values in a grid. Trond Kristiansen, 4.6.2009 on flight to Raleigh/Durham CO330
    """
    Nvars=len(varlist)
    t= (grdSTATION.endIndex-grdSTATION.startIndex)
    var_array_raw=np.zeros((t, int(len(grdSTATION.depth)), Nvars),dtype=np.float64)
    var_number=0
    
    grdSTATION.time=cdf.variables['time'][grdSTATION.startIndex:grdSTATION.endIndex]

    for var in varlist:
        print var_array_raw.shape
        var_array_raw[:,:,var_number] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex,:]
        
        """
        Slice the data in the netcdf file that you want to keep. In our case we store a time-series from
        a station (eta_rho, xi_rho), from top to bottom. The variables extracted at this station are
        defined in the "varlist" variable.
        """
        
        if log is True:
            print "\n---> Extracting time series of %s from station"%(var)
            print '---> Maximum %s %f'%(var, np.amax(var_array_raw[:,:,int(var_number)]))
            print '---> Minimum %s %f'%(var, np.amin(var_array_raw[:,:,int(var_number)]))      
        
        var_number=var_number+1
        
    grdSTATION.data=var_array_raw   
  