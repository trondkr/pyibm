import Nio
import numpy, datetime, types, math


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 10)
__modified__ = datetime.datetime(2008, 6, 10)
__version__  = "1.1"
__status__   = "Production"

def get_data(cdf_file, varlist, time, z_r, eta, xi, log):
    
    s_rho=cdf_file.variables["s_rho"].get_value()
    Nvars=len(varlist)
    var_array_raw=numpy.zeros((time, int(len(s_rho)), Nvars),float)
    var_number=0
    
    for var in varlist:
       
        
        tmp = cdf_file.variables[var].get_value()
        var_array_raw[:,:,var_number]=tmp[0:time,:,eta,xi]
        
        """
        Slice the data in the netcdf file that you want to keep. In our case we store a time-series from
        a station (eta_rho, xi_rho), from top to bottom. The variables extracted at this station are
        defined in the "varlist" variable.
        """
        
        if log is True:
            print "\n---> Extracting time series of %s from station (%s, %s)"%(var, xi, eta)
            print '---> Maximum %s %f'%(var, numpy.amax(var_array_raw[:,:,int(var_number)]))
            print '---> Maximum %s %f'%(var, numpy.amin(var_array_raw[:,:,int(var_number)]))      
        
        var_number=var_number+1
        
    return var_array_raw
 
def read_netcdf_file(filename, log):
    """
    Open the netCDF file and store the contents in arrays associated with variable names
    """
    try:
        cdf_file = Nio.open_file(filename,"r")
    except IOError:
        print 'Could not open file %s'%(filename)
        print 'Exception caught in: read_netcdf(filename, log)'
    if log is True:  
        print '\n---> Opened input file %s'%(filename)
    
    return cdf_file