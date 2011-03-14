
import netCDF4
import help
from netCDF4 import Dataset
import types, math
import numpy as np
import datetime as datetime
import os, sys
from memoize import memoize

""" Get self made modules"""
dirMine='/Users/trond/Projects/arcwarm/SODA/soda2roms'

if os.path.isdir(dirMine):
    sys.path.append(dirMine)

import grd
import IOstation
from initLarva import *

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 10)
__modified__ = datetime.datetime(2009, 12, 2)
__version__  = "1.2"
__status__   = "Production, 10.6.2008,2.12.2009, 24.06.2010"

def getStationData(cdf, varlist, grdSTATION, log, stationName,eventOption):
    """
    This routine reads a netCDF4 file created using IOstation.py in soda2roms. Such
    a file only contain information for one lat/long location and is therefore different
    from finding values in a grid. Trond Kristiansen, 4.6.2009 on flight to Raleigh/Durham CO330
    Edited in Bergen 24.06.2010. Edited again on 30.01.2011 to enable more variables from ESM2 runs.
    """
    if "useZooplanktonModel" in grdSTATION.OPTIONS:
        varlist.append("zoop")
    else:
        """Define the array as a scalar in the zooplankton variable in case not used."""
        grdSTATION.zooplankton=0.0

    Nvars=len(varlist)
    t= (grdSTATION.endIndex+1-grdSTATION.startIndex)
    var_numberXYZ=0
    var_numberXY =0
    """Create empty arrays based on number of 2D and 3D variables"""""
    for var in varlist:
        if var in ["taux","tauy","ssh"]:
            var_numberXY=var_numberXY+1
        else:
            var_numberXYZ=var_numberXYZ+1

    var_array_rawXY=np.zeros((t, var_numberXY),dtype=np.float64)
    var_array_rawXYZ=np.zeros((t, int(len(grdSTATION.depth)), var_numberXYZ),dtype=np.float64)

    var_numberXYZ=0
    var_numberXY =0
    """Fill the empty arrays with data and assign to grdSTATION object"""
    grdSTATION.time=cdf.variables['time'][grdSTATION.startIndex:grdSTATION.endIndex+1]

    for var in varlist:

        if var in ["taux","tauy","ssh"]:
            """2D variables"""
            var_array_rawXY[:,var_numberXY] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex+1]
            var_numberXY=var_numberXY+1

        elif var in ["zoop"]:
            """Special case where you calculate and add zooplankton based on other
            variables from the ESM files. This routine calls the Fortran module zooplankton.
            The final result is stored in a spedial array in grdSTATION.zooplankton as
            the array contains three dimensions (time,depth,prey item)."""
            help.showZooplanktonInfo()

            import zooplankton

            TEMP=np.squeeze(var_array_rawXYZ[:,:,0])
            NO3SM=np.squeeze(var_array_rawXYZ[:,:,4])
            NO3LG=np.squeeze(var_array_rawXYZ[:,:,5])
            NH4SM=np.squeeze(var_array_rawXYZ[:,:,2])
            NH4LG=np.squeeze(var_array_rawXYZ[:,:,3])
            TT=int(len(TEMP[:,0]))
            DD=int(len(TEMP[0,:]))
            II=int(len(prey_WGT))
            accZoopArray=np.zeros((TT,DD,II), dtype=np.float64)
            numberOfDays=30

            accZoopArray = zooplankton.zooplankton.calculatezoo(np.asarray(accZoopArray, order='Fortran'),
                                                                    NO3SM,NO3LG,
                                                                    NH4SM,NH4LG,
                                                                    TEMP,prey_D,
                                                                    prey_WGT,numberOfDays,TT,DD,II)
            grdSTATION.zooplankton=accZoopArray
        else:
            """3D variables"""
            var_array_rawXYZ[:,:,var_numberXYZ] = cdf.variables[var][grdSTATION.startIndex:grdSTATION.endIndex+1,:]
            var_numberXYZ=var_numberXYZ+1


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
            if var in ["taux","tauy","ssh"]:
                print '---> Maximum %s %3.16f'%(var, grdSTATION.dataXY[:,kXY].max())
                print '---> Minimum %s %3.16f'%(var, grdSTATION.dataXY[:,kXY].min())
                print '---> Mean %s %3.6f +-%3.16f'%(var, np.mean(grdSTATION.dataXY[:,kXY]), np.std(grdSTATION.dataXY[:,kXY]))
                kXY+=1
            elif var in ["zoop"]:
                for prey in range(len(grdSTATION.zooplankton[0,0,:])):
                    print '\t %i ---> Maximum %3.16f'%(prey, grdSTATION.zooplankton[:,0:index,prey].max())
                    print '\t    ---> Minimum %3.16f'%(grdSTATION.zooplankton[:,0:index,prey].min())
                    print '\t    ---> Mean %3.6f +-%3.16f \n'%(np.mean(grdSTATION.zooplankton[:,0:index,prey]), np.std(grdSTATION.zooplankton[:,0:index,prey]))
            else:
                print '---> Maximum %s %3.16f'%(var, grdSTATION.data[:,0:index,kXYZ].max())
                print '---> Minimum %s %3.16f'%(var, grdSTATION.data[:,0:index,kXYZ].min())
                print '---> Mean %s %3.6f +-%3.16f'%(var, np.mean(grdSTATION.data[:,0:index,kXYZ]), np.std(grdSTATION.data[:,0:index,kXYZ]))
                kXYZ+=1

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
        if (abs(ave.variables["temp"][:,latindex,lonindex,:]).any()< 10000):
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