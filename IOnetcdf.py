import netCDF4
from netCDF4 import Dataset
import types, math
import numpy as np
import datetime as dt
__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = dt.datetime(2008, 6, 10)
__modified__ = dt.datetime(2009, 12, 2)
__version__  = "1.2"
__status__   = "Production, 10.6.2008,2.12.2009, 24.06.2010"

def getStationData(cdf, varlist, grdSTATION, log, stationName):
    """
    This routine reads a netCDF4 file created using IOstation.py in soda2roms. Such
    a file only contain information for one lat/long location and is therefore different
    from finding values in a grid. Trond Kristiansen, 4.6.2009 on flight to Raleigh/Durham CO330
    Edited in Bergen 24.06.2010.
    """
    Nvars=len(varlist)
    t= (grdSTATION.endIndex+1-grdSTATION.startIndex)
    var_array_rawXY=np.zeros((t, Nvars),dtype=np.float64)
    var_array_rawXYZ=np.zeros((t, int(len(grdSTATION.depth)), Nvars),dtype=np.float64)
    
    var_numberXYZ=0
    var_numberXY =0
    
    grdSTATION.time=cdf.variables['time'][grdSTATION.startIndex:grdSTATION.endIndex+1]
   
    for var in varlist:
    
        if var in ["taux","tauy","ssh","chla"]:
            var_array_rawXY[:,var_numberXY] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex+1]
            var_numberXY=var_numberXY+1
        else:
            var_array_rawXYZ[:,:,var_numberXYZ] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex+1,:]
            var_numberXYZ=var_numberXYZ+1
                  
        """
        Slice the data in the netcdf file that you want to keep. In our case we store a time-series from
        a station (eta_rho, xi_rho), from top to bottom. The variables extracted at this station are
        defined in the "varlist" variable.
        """
        """Masked extreme values"""
   
    grdSTATION.data=var_array_rawXYZ
    grdSTATION.dataXY=var_array_rawXY
    grdSTATION.lat = float(cdf.variables['lat'][:])
   
    """Find and set the deepest depth for this station based on the depth range where you have
    valid data to use. Trond Kristiansen, 17.06.2010"""
    maxDepth=0; index=0
    for k in range(len(grdSTATION.data[0,:])):
        if ( abs(grdSTATION.data[0,k,0]) < 1000.):
            maxDepth = grdSTATION.depth[k]
            index=k
    maxDepth = round(maxDepth,0)
    grdSTATION.deepestDepthAllowed = maxDepth
    
    print "\nIOnetcdf.getData => Maximum depth for station %s has been set to: %s"%(stationName,grdSTATION.deepestDepthAllowed)
        
    if log is True:
        kXYZ=0;kXY=0
        for var in varlist:
            print "\n---> Extracted time series of %s from station %s"%(var,stationName)
            if var in ["taux","tauy","ssh","chla"]:
                print '---> Maximum %s %3.6f'%(var, grdSTATION.dataXY[:,kXY].max())
                print '---> Minimum %s %3.6f'%(var, grdSTATION.dataXY[:,kXY].min())
                print '---> Mean %s %3.6f +-%3.6f'%(var, np.mean(grdSTATION.dataXY[:,kXY]), np.std(grdSTATION.dataXY[:,kXY]))
                kXY+=1
            else:
                print '---> Maximum %s %3.6f'%(var, grdSTATION.data[:,0:index,kXYZ].max())
                print '---> Minimum %s %3.6f'%(var, grdSTATION.data[:,0:index,kXYZ].min())      
                print '---> Mean %s %3.6f +-%3.6f'%(var, np.mean(grdSTATION.data[:,0:index,kXYZ]), np.std(grdSTATION.data[:,0:index,kXYZ]))
                kXYZ+=1
            
def getData(julian,julianIndex,julianFileA,julianFileB,dz1,dz2,depthIndex1,depthIndex2,grdSTATION,event):
    """Calculate weights to use on input data from file"""
    dwB = abs(julian) - abs(julianFileA)
    dwA = abs(julianFileB) - abs(julian)
   
    if event != "ESM RUN":
        """Interpolate the values of temp, salt, u and v velocity in time to current julian date"""
        Tdata=((grdSTATION.data[julianIndex,depthIndex1,0])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,0])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,0])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,0])*
            (dwB/(dwA+dwB)))*dz2
        Sdata=((grdSTATION.data[julianIndex,depthIndex1,1])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,1])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,1])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,1])*
            (dwB/(dwA+dwB)))*dz2
        Udata=((grdSTATION.data[julianIndex,depthIndex1,2])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,2])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,2])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,2])*
            (dwB/(dwA+dwB)))*dz2
        Vdata=((grdSTATION.data[julianIndex,depthIndex1,3])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,3])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,3])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,3])*
            (dwB/(dwA+dwB)))*dz2
        TauX=((grdSTATION.dataXY[julianIndex,0])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,0])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.dataXY[julianIndex,0])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,0])*
            (dwB/(dwA+dwB)))*dz2
        TauY=((grdSTATION.dataXY[julianIndex,1])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,1])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.dataXY[julianIndex,1])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,1])*
            (dwB/(dwA+dwB)))*dz2
        windX, windY = convertStressToWind(TauX,TauY)
        #print '---> Maximum %s %f'%("windX", windX.max())
        #print '---> Minimum %s %f'%("windY", windY.min())
    
        return Tdata,Sdata,Udata,Vdata, windX, windY

    else:
        """Interpolate the values of temp, salt, u and v velocity in time to current julian date.
        Note that the indices for 2D and 3D values are different and must be set accordingly to
        the order the 2D and 3D variables are stored into arrays. The order is defined by the order they appear
        in vars list defined in init funtion."""
        Tdata=((grdSTATION.data[julianIndex,depthIndex1,0])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,0])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,0])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,0])*
            (dwB/(dwA+dwB)))*dz2
        NSMdata=((grdSTATION.data[julianIndex,depthIndex1,1])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,1])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,1])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,1])*
            (dwB/(dwA+dwB)))*dz2
        NLGdata=((grdSTATION.data[julianIndex,depthIndex1,2])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex1,2])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.data[julianIndex,depthIndex2,2])*
            (dwA/(dwA+dwB))+(grdSTATION.data[julianIndex+1,depthIndex2,2])*
            (dwB/(dwA+dwB)))*dz2
        CHLAdata=((grdSTATION.dataXY[julianIndex,0])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,0])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.dataXY[julianIndex,0])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,0])*
            (dwB/(dwA+dwB)))*dz2
        TauX=((grdSTATION.dataXY[julianIndex,1])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,1])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.dataXY[julianIndex,1])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,1])*
            (dwB/(dwA+dwB)))*dz2
        TauY=((grdSTATION.dataXY[julianIndex,2])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,2])*
            (dwB/(dwA+dwB)))*dz1 +((grdSTATION.dataXY[julianIndex,2])*
            (dwA/(dwA+dwB))+(grdSTATION.dataXY[julianIndex+1,2])*
            (dwB/(dwA+dwB)))*dz2
        windX, windY = convertStressToWind(TauX,TauY)
      #  print '---> Maximum %s %f'%("windX", windX.max())
      #  print '---> Minimum %s %f'%("windY", windY.min())
        
        return Tdata,NSMdata,NLGdata,CHLAdata,windX,windY


def convertStressToWind(TauX,TauY):
    """Convert surface stress to wind velocity for use in turbulence calculations of
    encounter rate between larvae and prey. We assume a drag coefficient
    with the shape Cd=0.44+0.063U (where U is wind at 10 m).
    In SODA TauX and TauY are in N/m**2. This is equivavlent to
    N=kgm/s**2 and therefore kg/ms**2. The unit for air density is kg/m**3 so
    after doing the quadratic estimate of U we get wind velocity with units m/s.
    
    Update: changed drag coefficient top constant = 1.2e-3 and wind density to 1.3. This sfunction also
    disregards direction of wind as we are only concerned with the scalar value."""
    windX=abs(np.sqrt(abs(TauX)/(1.3*1.2e-3)))
    windY=abs(np.sqrt(abs(TauY)/(1.3*1.2e-3)))
   
    return windX, windY
    
def getAverageValues(inFile,lon,lat):
    """inFile is the climatology file (e.g. averageSODA1961-1990.nc), while lat and lon is the
    station longitude/latitude"""
    numberOfPoints=4
    grdMODEL = grd.grdClass(inFile,"AVERAGE")
    
    gridIndexes, dis = IOstation.getStationIndices(grdMODEL,lon,lat,'AVERAGE',numberOfPoints)
    ave = Dataset(inFile,'r')
    #for t in range(12):
    #    test=ave.variables["temp"][1,:,:,t]
    #    test=np.ma.masked_less(test,-1000)
    #    contourf(np.squeeze(test),100)
    #    show()
    aveTemp=np.zeros((12),np.float64)
    #aveSalt=np.zeros((12),np.float64)
    #aveUvel=np.zeros((12),np.float64)
    #aveVvel=np.zeros((12),np.float64)
    
    validIndex=[]; validDis=[]
    for i in range(numberOfPoints):
        latindex=int(gridIndexes[i][0])
        lonindex=int(gridIndexes[i][1])
        
        """Test to see if station contains anything but missing values. If it does, we assume
        this holds for salinity, u, and v as well. We also assume all values less than 10000"""
        if any(abs(ave.variables["temp"][:,latindex,lonindex,:])< 10000):
            validIndex.append(i)
            validDis.append(dis[i])
  
    if not validIndex:
        print 'No valid data found for position'
        exit()
    
    for i in validIndex:
        for month in range(12):
            wgt=float(validDis[i])/sum(validDis)
            latindex=int(gridIndexes[i][0])
            lonindex=int(gridIndexes[i][1])
            
            """The values at a station is calculated by interpolating from the
            numberOfPoints around the station uysing weights (wgt)
            """
            aveTemp[month] = aveTemp[month]  + np.mean((ave.variables["temp"][0:2,latindex,lonindex,month])*wgt)
            #aveSalt[month] = aveSalt[month]  + np.mean(ave.variables["salt"][0:2,latindex,lonindex,month])*wgt
            #aveUvel[month] = aveUvel[month]  + np.mean(ave.variables["u"][0:2,latindex,lonindex,month])*wgt
            #aveVvel[month] = aveVvel[month]  + np.mean(ave.variables["v"][0:2,latindex,lonindex,month])*wgt
    print "Finished storing average values for station into array"   
    ave.close()
    
    return aveTemp

def getSeaWifs(seawifsFileName):
    """Note that the oder of the station data read from seawifs file needs to
    correspond to the order of stations calulcated and defined in main funtion!!!"""
    file = Dataset(seawifsFileName)
    print "SEAWIFS data -------------------------------------------"
    print "Extracting seawifs data from file %s"%(seawifsFileName)
    print "%s"%(file.stations)
    seawifsArray = file.variables['Chlorophyll-a']
   
    return seawifsArray