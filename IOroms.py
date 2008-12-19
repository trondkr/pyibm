import Nio
import numpy, math, datetime

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 13)
__modified__ = datetime.datetime(2008, 6, 13)
__version__  = "1.0"
__status__   = "Development"


def sigma2meters(cdf_file, time, log):
    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[eta_rho, xi_rho, s]
    This calculation takes time and should be matrix manipulated to speed it up.
    Trond Kristiansen, 20.01.2008
    """
    
    Cs_r=cdf_file.variables["Cs_r"].get_value() 
    Cs_w=cdf_file.variables["Cs_w"].get_value() 
    h   =cdf_file.variables["h"].get_value()  
    hc  =cdf_file.variables["hc"].get_value()  
    zeta=cdf_file.variables["zeta"].get_value() 
    s_rho=cdf_file.variables["s_rho"].get_value()   
    
    dim=cdf_file.variables["zeta"].shape
 
    eta_rho=dim[1]
    xi_rho=dim[2]
  
    """
    Crearting a new array (z_r) that store the sigma to meters matrix
    """
    z_r=numpy.zeros((time,len(s_rho),eta_rho,xi_rho),float)
    zeta_edit=numpy.zeros((eta_rho,xi_rho),float)
   
    for t in range(time):
     for k in range(len(s_rho)):
        zeta_edit=zeta[t,:,:]
        z_r[t,k,:,:] = numpy.multiply(zeta_edit,(1+s_rho[k]))+hc*(s_rho[k])+(numpy.subtract(h,hc))*Cs_r[k]
         
        
    if log is True:
        
        print 'Minimum depth %s'%(numpy.amax(z_r))
        print 'Maximum depth %s'%(numpy.amin(z_r))
        print 'Number of sigma-layers in file %s'%(len(s_rho))
        print 'Saved data for %s timesteps'%(time)
        
    return z_r

