import os, sys, string
from shutil import move
import netCDF4
from netCDF4 import num2date, Dataset
import types
import numpy as np
import IOwrite
import datetime as datetime
from pylab import *
"""
Import modules created as part of project
"""

""" Get self made modules"""
dirMine='/Users/trond/Projects/arcwarm/SODA/soda2roms'

if os.path.isdir(dirMine):
    sys.path.append(dirMine)

import grd
import IOverticalGrid
import IOtime
#import IOlight
#import date
import IOnetcdf
import predation
from initLarva import *
"""Import f2py Fortran modules"""
import calclight
import perception
import IOstation

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 9)
__modified__ = datetime.datetime(2010, 2, 10)
__version__  = "1.1.5"
__status__   = "Production"
__revisions__= "8.6.09, 2.12.09, 19.01.10, 26.01.10, 10.02.2010"

    
def infoOnResolution(grdSTATION):
  
    print 'Defined resolution of simulation:'
    print 'Temporal : %s hours per time-steps (%s hours)'%(deltaH, grdSTATION.refDate +
                                                             datetime.timedelta(seconds=deltaH*3600.)
                                                             - grdSTATION.refDate)
    print 'Vertical : %s cm (rounds to %i decimalpoints)'%(deltaZ*100.,lastDecimal)
    print 'The simulations are run with %i prey densities'%(Nprey)
    print '-----------------------------------------------------------------\n'

def coldWarmEvent(stationName,event):
    #['North Sea','Iceland','East Greenland', 'Lofoten', 'Georges Bank']
    coldYears=[1963,1983,1962,1977,1965]
    warmYears=[1999,1960,1963,1990,1999]
    
    if event=='COLD':
        if stationName=='North Sea':
            startDate = datetime.datetime(coldYears[0],1,15,0,0,0)
            endDate   = datetime.datetime(coldYears[0],12,31,0,0,0)
        if stationName=='Iceland':
            startDate = datetime.datetime(coldYears[1],1,15,0,0,0)
            endDate   = datetime.datetime(coldYears[1],12,31,0,0,0)
        if stationName=='East Greenland':
            startDate = datetime.datetime(coldYears[2],1,15,0,0,0)
            endDate   = datetime.datetime(coldYears[2],12,31,0,0,0)
        if stationName=='Lofoten':
            startDate = datetime.datetime(coldYears[3],1,15,0,0,0)
            endDate   = datetime.datetime(coldYears[3],12,31,0,0,0)
        if stationName=='Georges Bank':
            startDate = datetime.datetime(coldYears[4],1,15,0,0,0)
            endDate   = datetime.datetime(coldYears[4],12,31,0,0,0)
    if event=='WARM':
        if stationName=='North Sea':
            startDate = datetime.datetime(warmYears[0],1,15,0,0,0)
            endDate   = datetime.datetime(warmYears[0],12,31,0,0,0)
        if stationName=='Iceland':
            startDate = datetime.datetime(warmYears[1],1,15,0,0,0)
            endDate   = datetime.datetime(warmYears[1],12,31,0,0,0)
        if stationName=='East Greenland':
            startDate = datetime.datetime(warmYears[2],1,15,0,0,0)
            endDate   = datetime.datetime(warmYears[2],12,31,0,0,0)
        if stationName=='Lofoten':
            startDate = datetime.datetime(warmYears[3],1,15,0,0,0)
            endDate   = datetime.datetime(warmYears[3],12,31,0,0,0)
        if stationName=='Georges Bank':
            startDate = datetime.datetime(warmYears[4],1,15,0,0,0)
            endDate   = datetime.datetime(warmYears[4],12,31,0,0,0)
    print 'Special climate run with %s event define: start %s for station %s'%(event,startDate,stationName)
    return startDate, endDate

def init(station,stationName,event):
    """
    init: initializes the reading of files and definition of global variables
    """
    log         = True
    clim        = False
    Finish      = False
    useAverageFile = True
    fileNameIn  = station

    """Option for using seawifs data as prey abundance"""
    seawifs=True
    seawifsFileName="/Users/trond/Projects/seawifs/chlo-stations.nc"
    if seawifs==True:
        seawifsArray=getSeaWifs(seawifsFileName)
    else:
        seawifsArray=None
        
    """Open the netcdf file if it existst."""
    cdf = Dataset(fileNameIn)
  
    """Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function."""
    grdSTATION = grd.grdClass(fileNameIn,"STATION")
    
    if useAverageFile==True:
        inFile="/Users/trond/Projects/arcwarm/SODA/soda2average/clim/averageSODA1961-1990.nc"
        aveTemp = getAverageValues(inFile,grdSTATION.lon,grdSTATION.lat)
        grdSTATION.aveT=aveTemp
    
    varlist=['temp','salt','u','v','taux','tauy'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']

    startDate = datetime.datetime(1970,4,1,0,0,0)
    endDate   = datetime.datetime(1970,4,6,0,0,0)
    if event == 'COLD' or event=='WARM':
        startDate, endDate = coldWarmEvent(stationName,event)
        
    """Get the time information and find the indices for start and stop data to extract relative to
    the time period wanted. Takes input data:
    startDate="DD/MM/YYYY" and endDate="DD/MM/YYYY" or if none given, finds all date"""
    getTimeIndices(cdf,grdSTATION,startDate,endDate,clim)
    
    """Get the total number of hours so that we can use this number to create arrays."""
    delta=(endDate - startDate)
    grdSTATION.delta=int((delta.days + 1./(24./deltaH))*(24./deltaH))
    grdSTATION.endDate=endDate
    grdSTATION.startDate=startDate
        
    if grdSTATION.delta/24*deltaH < NDaysAlive:
        print '\nWARNING: You are running for time-period less than days you want cohort to stay alive'
        print 'Try to change NDaysAlive from %s to %s'%(NDaysAlive,int(grdSTATION.delta/24*deltaH))
        print 'You can change this in initLarva.py'
        sys.exit()
    
    """Extract the variables at the given station and store in the grdSTATION object"""
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,log,clim,stationName)
    """Store the maximum and minimum temperature values for the time period of interest, and
    use these to find the prey density as a function of temperature: function calculateGrowth"""
    
    grdSTATION.minT = grdSTATION.data[:,:,0].min()
    grdSTATION.maxT = grdSTATION.data[:,:,0].max()
    
    """Open output file:"""
    outputFile="IBM_"+str(startDate.year)+"_"+os.path.basename(station)
    if not os.path.exists("results"): os.mkdir("results")
    if os.path.exists(outputFile): os.remove(outputFile)
    if os.path.exists("results/"+outputFile): os.remove("results/"+outputFile)
    
    listOfReleaseDates=[]
    """Calculate release dates for individual cohorts based on days since start date of simulations."""
    listOfReleaseDates.append(startDate)
  
    for i in range(NReleaseDatesInit):
        date=startDate+datetime.timedelta(daysBetweenReleases*(i+1))
        if endDate >= date + datetime.timedelta(days=NDaysAlive):
            listOfReleaseDates.append(date)
        
    NReleaseDates=len(listOfReleaseDates)
    grdSTATION.listOfReleaseDates=listOfReleaseDates
    
    print "\nThis simulation will release a total of %s cohorts on the following dates:"%(NReleaseDates)
    for h in range(NReleaseDates):
        print "--> %s"%(listOfReleaseDates[h])
    print "\n"
    
    print "IBM will run from date %s to %s\n"%(grdSTATION.listOfReleaseDates[0] + datetime.timedelta(days=0),
                                              grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive))
    print "-----------------------------------------------------------------"                       
    
    grdSTATION.seawifs=seawifs
    grdSTATION.stationNumber=0
    
    infoOnResolution(grdSTATION)
      
    return grdSTATION, clim, outputFile, NReleaseDates, Finish, seawifsArray

def findActiveCohorts(isDead,isReleased):
    co=0; co2=0
    for check in range(len(isDead)):
        if isDead[check]==True: co+=1
   
    for check in range(len(isReleased)):
        if isReleased[check]==True: co2+=1
    activeCohorts=co2-co
    minNumberOfActiveCohort = co
   # print "There are currently %s active cohorts in the simulation with min cohort %s"%(activeCohorts,minNumberOfActiveCohort
   
    return activeCohorts, minNumberOfActiveCohort              
                  
def checkReleased(releaseDate,julian,grdSTATION):
  
    if grdSTATION.firstBatch==True:
        release=True
        isReleased=False
        grdSTATION.firstBatch=False
        releaseDate=grdSTATION.refDate + datetime.timedelta(seconds=julian)
    else:
        if releaseDate == grdSTATION.refDate + datetime.timedelta(seconds=julian):
            release = True
            isReleased= False
        if releaseDate > grdSTATION.refDate + datetime.timedelta(seconds=julian):
            release = False
            isReleased= False
              
        if releaseDate < grdSTATION.refDate + datetime.timedelta(seconds=julian):
            release = False
            isReleased= True
   
    #print "Will release:",release, " Has been released",isReleased, " at date:",releaseDate," current date:", grdSTATION.refDate + datetime.timedelta(seconds=julian)
   # print  "current date:", grdSTATION.refDate + datetime.timedelta(seconds=julian)
    return release, isReleased, releaseDate
    
    
def isAlive(julian,larvaAge,cohort,ind,t,prey,startAndStopIndex,isDead,grdSTATION):
 
    if isDead[cohort]==False:
       # print "Larva is %s days old and belongs to cohort %s %s"%(larvaAge/24.,cohort, NDaysAlive)
        
        if larvaAge/24. > int(NDaysAlive):
            startAndStopIndex[cohort,:,1]=t-24
            isAliveBool=False
            isDead[cohort]=True
        else:
            isAliveBool=True
    else:
        isAliveBool=False
    grdSTATION.startAndStopIndex=startAndStopIndex
   
    return isAliveBool, grdSTATION
    
def getSeaWifs(seawifsFileName):
    """Note that the oder of the station data read from seawifs file needs to
    correspond to the order of stations calulcated and defined in main funtion!!!"""
    file = Dataset(seawifsFileName)
    print "SEAWIFS data -------------------------------------------"
    print "Extracting seawifs data from file %s"%(seawifsFileName)
    print "%s"%(file.stations)
    seawifsArray = file.variables['Chlorophyll-a']
   
    return seawifsArray
    
def getDepthIndex(grdSTATION,depth):
    depthFOUND=False; d=0
    while d < len(grdSTATION.depth) and depthFOUND==False:
        
        if abs(grdSTATION.depth[d]) < abs(depth) <= abs(grdSTATION.depth[d+1]) and depthFOUND==False: 
            dz1=1.0-abs((abs(grdSTATION.depth[d])-abs(depth))/(abs(grdSTATION.depth[d])-abs(grdSTATION.depth[d+1])))
            dz2=1.0-abs((abs(grdSTATION.depth[d+1])-abs(depth))/(abs(grdSTATION.depth[d])-abs(grdSTATION.depth[d+1])))
            depthIndex1=d; depthIndex2=d+1
            depthFOUND=True
            #print 'depth 1:',grdSTATION.depth[d], depth, grdSTATION.depth[d+1],depthIndex1, depthIndex2
            
        elif abs(grdSTATION.depth[0]) > abs(depth) and depthFOUND==False:
            depthIndex1=0; depthIndex2=0; dz1=0.5; dz2=0.5
            depthFOUND=True
            #print 'depth 2:',grdSTATION.depth[d], depth,depthIndex1, depthIndex2
             
        elif abs(grdSTATION.depth[-1]) < abs(depth) and depthFOUND==False:
            depthIndex1=-1; depthIndex2=-1; dz1=0.5; dz2=0.5
            depthFOUND=True
            #print 'depth 3:',grdSTATION.depth[-1], depth,depthIndex1, depthIndex2
        elif d==len(grdSTATION.depth)-1 and depthFOUND==False:
            print "No valid index for depth (%sm) found"%(depth)
           
        d+=1
    return depthIndex1, depthIndex2, dz1, dz2

def initBehavior(depth,maxHourlyMove):
    oldDepth=depth  #The previous depth for a given larvae for a time step inside optimal loop
    optDepth=depth  # The optimal depth for a given larvae for a time step 
    h_start, h_stop = getMaxMinDepths(oldDepth,maxHourlyMove)
    NOptDepths=int((abs(oldDepth-h_start) + abs(h_stop-oldDepth))/float(deltaZ))
    fitness=-9999.9
  
    return oldDepth,optDepth,h_start,h_stop,NOptDepths,fitness

def getBehavior(stomachFullness,F,m,length,oldFitness,depth,optDepth):
    """Rule 4 of Behavioral Ecology paper - Kristiansen et al. 2009"""
    T=min(1.0,0.3+1000.0*(1+(length)*np.exp(length))**(-1))

    #print "behavior ", T, length, 0.3+1000.0*(1+(length)*np.exp(length))**(-1), length
    if stomachFullness > T:
        beta   = 7.0
        vector = ((stomachFullness-T)/(1.0-T))**beta
    else: vector=0.0
    
    fitness=(((1-vector)*F) - vector*m)
   
    #print "Vector %s optdepth %s vs old depth %s : new fit %s and old %s "%(vector,optDepth,depth,fitness,oldFitness) 
    if fitness > oldFitness:
     #   print "Change: Vector %s new depth %s vs old depth %s : new fit %s and old %s "%(vector,optDepth,depth,fitness,oldFitness) 
        optDepth=depth
        oldFitness=fitness
    
    return optDepth,oldFitness

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
    
def getData(julian,julianIndex,julianFileA,julianFileB,dz1,dz2,depthIndex1,depthIndex2,grdSTATION):
    """Calculate weights to use on input data from file"""
    dwB = abs(julian) - abs(julianFileA)
    dwA = abs(julianFileB) - abs(julian)
    #print dwA/(dwA+dwB),dwB/(dwA+dwB),julian, julianFileA,julianFileB
   
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

def calculateLength(Larval_wgt, Larval_m):
    oldLarval_m = Larval_m
    Larval_m = np.exp(2.296 + 0.277*np.log(Larval_wgt) - 0.005128*np.log(Larval_wgt)**2)
    if Larval_m < oldLarval_m:
        Larval_m = oldLarval_m 
        #print 'Kept the old length', oldLarval_m, Larval_m
    return Larval_m

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
    
def calculateGrowth(julian,Eb,deltaH,depth,hour,grdSTATION,dt,Larval_wgt, Larval_mm,
                    julianIndex,julianFileA,julianFileB,
                    R,enc,hand,ing,pca,Spre,prey):
    
    """==>INDEX calculations==============================="""
    depthIndex1, depthIndex2, dz1, dz2 = getDepthIndex(grdSTATION,depth)
    Tdata,Sdata,Udata,Vdata,windX,windY = getData(julian,julianIndex,julianFileA,
                                      julianFileB,dz1,dz2,depthIndex1,
                                      depthIndex2,grdSTATION)
                        
    """Mouthsize of larvae. The larvae can only capture calanus of sizes less
    than the larval mouthsize. Folkvord et al."""
    m = np.exp (-3.27+1.818*np.log(Larval_mm)-0.1219*(np.log(Larval_mm))**2.)
    """Calculate metabolism"""
    meta = dt*2.38e-7*np.exp(0.088*7)*((Larval_wgt*mg2ug)**(0.9)*0.001)*deltaH
    """Increase metabolism during active hours (light above threshold) Lough et  al. 2005"""
    if Eb > 0.001:
        if Larval_mm > 5.5: meta = (2.5*meta)
        else: meta = (1.4*meta)
    """Calculate assimilation efficiency"""
    assi=0.8*(1-0.400*exp(-0.002*(Larval_wgt*mg2ug-50.0)))
    """alculate growth rate"""
    GR = max(0.0, sec2day*np.log(0.01*(1.08 + 1.79*Tdata - 0.074*Tdata*np.log(Larval_wgt)
                              - 0.0965*Tdata*np.log(Larval_wgt)**2
                              + 0.0112*Tdata*np.log(Larval_wgt)**3) + 1)*dt*deltaH)
    """Growth rate as percent"""
    GR_gram = (np.exp(GR) - 1)*Larval_wgt
         
  
    """FOOD == Calculate seasonal prey density based on temperature and phytoplankton"""
    currentDate = grdSTATION.refDate + datetime.timedelta(seconds=julian)
    Tanomaly = (Tdata - grdSTATION.aveT[int(currentDate.month-1)]) / (grdSTATION.aveT.max() - grdSTATION.aveT.min())
    
    if grdSTATION.relativeSeawifsValue > 0.01:
        P = Tanomaly + grdSTATION.relativeSeawifsValue
    else:
        P = grdSTATION.relativeSeawifsValue
    
    P = max(0.0,P)
    
    """VISUAL == PERCEPTION of PREY calculations==============================="""          
    Em = (Larval_mm**2.0)/(contrast*0.1*0.2*0.75) #Size-specific sensitivity of the visual system (Fiksen,02)
    
    #for j in range(6):
    j=1
    IER=0; R[j]=0.0 #R[j]=np.sqrt(Em*contrast*(calanus_Area[j])*(Eb/(Ke_larvae+Eb)))
    """All input to getr is either in m (or per m), or in mm (or per mm)"""
    R[j], IER = perception.perception.getr(R[j],beamAttCoeff/m2mm,contrast,calanus_Area[j],Em,Ke_larvae,Eb, IER)
   
    """If you don't have fortran compiler use python routine for visual range"""
    #visual = np.sqrt(Em*contrast*(Ap_calanus*m2mm)*(Eb/(Ke_larvae+Eb)))
    #R[j] = (IOlight.getPerceptionDistance(Em,attCoeff,Ke_larvae,Ap_calanus,Eb))*m2mm
    
    """Calculate turbulence based on wind stress"""
    epsilon=(5.82*1.E-9*((np.sqrt(windX**2 + windY**2)))**3.)/(depth+0.1)
    omega = 1.9*(epsilon*R[j]*mm2m)**0.667
    omega =  omega * m2mm # From m/s to mm/s
    
    """Calculate handling time, encounter rate, and probability of capture"""
    hand[j] = 0.264*10**(7.0151*(calanus_L1[j]/Larval_mm)) # Walton 1992
    enc[j] = ((0.667*np.pi*(R[j]**3.)*f + np.pi*(R[j]**2.)*np.sqrt(calanus_L1[j]**2.+ 2.*omega**2.)*f*tau) * (calanus_D[j]*((prey+1)*MultiplyPrey*P))* ltr2mm3)
    
    
    pca[j] = enc[j]*max(0.0,min(1.0,-16.7*(calanus_L1[j]/Larval_mm) + 3.0/2.0))
    
    """Calculate ingestion rate"""
    ing[j] = (dt*pca[j]*calanus_W[j]*micro2m / (1 + hand[j]))*deltaH
   
    """Calculate stomach fullness"""
    stomachFullness =  (min(gut_size*Larval_wgt,Spre + sum(ing)))/(Larval_wgt*gut_size)
  
    return ing,GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness, MultiplyPrey*(prey+1)*P
   
def getMaxMinDepths(oldDepth,maxHourlyMove):
    """The value you give deepstDepthAllowed will be the absolute deepest
    depth the larvae will move to. Useful to set this depth equal to bottom."""
    deepstDepthAllowed=55.
    if (oldDepth>maxHourlyMove and oldDepth+maxHourlyMove < deepstDepthAllowed):
        h_start = oldDepth - maxHourlyMove
        h_end   = oldDepth + maxHourlyMove

    elif (oldDepth <= maxHourlyMove):
        h_start=0.0
        h_end=min(deepstDepthAllowed,oldDepth + maxHourlyMove)
    elif (oldDepth + maxHourlyMove >= deepstDepthAllowed):
        h_start=max(1,oldDepth-maxHourlyMove)
        h_end=deepstDepthAllowed
    else:
        print 'STOP: Depth error (ibm.py:getMaxMinDepths)'
    
    return h_start, h_end


def calculateAttenuationCoeff(Chla):
    """See Fiksen et al. 2002 for details (relationship taken from Riley 1956)
    units of Chla is mg/m3"""
    k0 = 0.1
    return k0 + 0.054*Chla**(2./3.) + 0.0088*Chla
    
def getTimeIndices(cdf,grdSTATION,startDate=None,endDate=None,clim=None):
    """Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    refDate    = datetime.datetime(1948,1,1,0,0,0)
    grdSTATION.refDate = refDate
    
    if clim is False:
        """Get start and stop time stamps from file. The time stamp from SODA is
        in days so conver to seconds (since 1948,1,1)"""
        startJDFile =int(grdSTATION.time[0])*86400
        endJDFile   =int(grdSTATION.time[-1])*86400
        
        """Figure out the date from the JD number. I have consistently used JD numbers
        relative to 01/01/1948 and you therefore have to add that jdref number
        to calculate the correct date using the date module."""
     
        startDateFile = refDate + datetime.timedelta(seconds=startJDFile)
        endDateFile   = refDate + datetime.timedelta(seconds=endJDFile)
       
        print "\nStation contains data for the time period:"
        print "%s to %s"%(startDateFile, endDateFile)
        
        """Now find the Juian dates of the time period you have asked for and see
        if it exists in the file. If so, the find the indices that correspond to the
        time period.
        """
        JDstart = startDate - refDate
        JDstart = JDstart.days*86400 + JDstart.seconds
        
        JDend   = endDate - refDate
        JDend   = JDend.days*86400 + JDend.seconds
    
        if JDstart < startJDFile:
            print "Start time preceeds the earliest time stamp found in file"
            print "Required in File: %s"%(refDate + datetime.timedelta(seconds=JDstart))
            print "Actually in File: %is"%(refDate + datetime.timedelta(seconds=startJDFile))
            exit()
            
        if JDend > endJDFile:
            print "End time exceeds the last time stamp found in file"
            print "Required in File: %s"%(refDate + datetime.timedelta(seconds=JDend))
            print "Actually in File: %is"%(refDate + datetime.timedelta(seconds=endJDFile))
            exit()
            
        if JDend < endJDFile and JDstart > startJDFile:
            print "Time period to extract was found within the time period available in the file..."
            print "--> %s - %s"%(startDate,endDate)
            
            for i in range(grdSTATION.time.shape[0]):
               
                if grdSTATION.time[i]*86400  < JDstart:
                    FOUND=False
                    continue
                elif grdSTATION.time[i]*86400  >= JDstart and FOUND is False:
                    print "\nFound first time index that fits start point:"
                    print "%s at index %i => %s"%(refDate + datetime.timedelta(seconds=grdSTATION.time[i-1]*86400),i-1,JDstart)
                    grdSTATION.startIndex=i-1
                    grdSTATION.start_date = refDate + datetime.timedelta(seconds=grdSTATION.time[i-1]*86400)
                    FOUND=True
                    
            for i in range(grdSTATION.time.shape[0]):
                if grdSTATION.time[i]*86400 < JDend:
                    FOUND=False
                    continue
                elif grdSTATION.time[i]*86400 >= JDend and FOUND is False:
                    print "Found first time index that fits end point:"
                    print "%s at index %i => %s\n"%(refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400),i,JDend)
                    grdSTATION.endIndex=i+1
                    grdSTATION.end_date = refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400)
                    FOUND=True
        print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        
        grdSTATION.JDend=JDend
        grdSTATION.JDstart=JDstart
        grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        
    else:
        print "Using climatological time as clim is set to True"
        d = string.split(startDate,'/')
        start_date = date.Date()
        start_date.day=int(d[0])
        start_date.month=int(d[1])
        start_date.year=int(d[2])
        jdstart=start_date.GetYearDay()
        
        d = string.split(endDate,'/')
        end_date = date.Date()
        end_date.day=int(d[0])
        end_date.month=int(d[1])
        end_date.year=int(d[2])
        jdend=end_date.GetYearDay()
        
        grdSTATION.time=[]
        grdSTATION.time=(cdf.variables["clim_time"][:])*5
        FOUND = False
        
        for i in range(grdSTATION.time.shape[0]):
           
            if grdSTATION.time[i]  < jdstart:
                continue
            elif grdSTATION.time[i]  >= jdstart and FOUND is False:
                print "\nFound first time index that fits start point:"
                print "%s at index %i => %s"%(grdSTATION.time[i],i,jdstart)
                grdSTATION.startIndex=i-1
                grdSTATION.start_date = date.Date()
                grdSTATION.start_date =date.DateFromJDNumber(jdstart)
                FOUND=True
        for i in range(grdSTATION.time.shape[0]):
            if grdSTATION.time[i] < jdend:
                FOUND=False
                continue
            elif grdSTATION.time[i]  >= jdend and FOUND is False:
                print "Found first time index that fits end point:"
                print "%s at index %i => %s\n"%(grdSTATION.time[i],i,jdend)
                grdSTATION.endIndex=i+1
                grdSTATION.end_date = date.Date()
                grdSTATION.end_date =date.DateFromJDNumber(jdend)
                FOUND=True
                
        grdSTATION.jdstart=jdstart
        grdSTATION.jdend=jdend
        print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])
        
def updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION,clim):
    """This routine makes sure that we are reading at the correct time indices
    in the input data. If we have moved forward in time enough to move past the current time window of two
    time stamps, then update the indices and give the new indices"""
   
    if julian >=julianFileB and julian < grdSTATION.JDend:
        julianIndex=julianIndex+1
        julianFileA=grdSTATION.time[julianIndex]*86400
        julianFileB=grdSTATION.time[julianIndex+1]*86400
        Finish = False
    elif julian < julianFileB and julian < julianFileA:
        """Now, if you run with more than one individual you have to reset the
        fileIndices for where in the input files you interpolate from between each individual.
        Therefore, here we search backwards to find the correct start and stop indices for the new
        individual."""
        while julian < julianFileB and julian < julianFileA:
            julianIndex=julianIndex-1
            julianFileA=grdSTATION.time[julianIndex]*86400
            julianFileB=grdSTATION.time[julianIndex+1]*86400
            Finish = False
            
    elif julian > grdSTATION.JDend:
        Finish = True
    else:
        Finish = False
    if clim is False:
        julianFileA=grdSTATION.time[julianIndex]*86400
        julianFileB=grdSTATION.time[julianIndex+1]*86400
    else:
        julianFileA=grdSTATION.time[julianIndex]
        julianFileB=grdSTATION.time[julianIndex+1]
    
    return julianFileA, julianFileB, julianIndex, Finish

def showProgress(t,Ntime,message):
    p=int( ((t*1.0)/(1.0*(Ntime-24)))*100.)
    progress.render(p,message)

def ibm(station,stationName,stationNumber,event):
    
    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION, clim, outputFile, Ncohorts, Finish, seawifsArray = init(station,stationName,event)
    
    """Start time in days, hours, and seconds since 1948,1,1"""
    time=grdSTATION.startDate - grdSTATION.refDate 
   
    """ julian is the time used to control the larvae, while
    julianFileA and julianFileB are the time stamps of inbetween
    the larvae recides relative to time from file. Temperature is
    interpolated to date julian from julianFileA and julianFileB"""
    julian=time.days*86400 + time.seconds
    
    julianFileA=grdSTATION.time[0]*86400
    julianFileB=grdSTATION.time[1]*86400
    Ntime=int(grdSTATION.delta) # Number of time-steps to use in arrays
        
    if clim is True:
        julianFileA=grdSTATION.time[0]
        julianFileB=grdSTATION.time[1]
    julianIndex=0
    
    """Create array for each individual larvae"""
   
    W   =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    larvaPsur   =(np.ones((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))
    larvaDepth   =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    W_AF=(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    S   =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    Age =(np.ones((Ncohorts,Nlarva,Nprey),dtype=np.float64))*missingValue
    startAndStopIndex =(np.ones((Ncohorts,Nlarva,2),dtype=np.float64))*missingValue
    L   =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
   
    ingrate  =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    Psurvive =(np.ones((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    SGR      =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    larvaTdata   =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
    larvaNauplii  =(np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64))*missingValue
   
    R=np.zeros(13,dtype=np.float64)
    enc=np.zeros(13,dtype=np.float64)
    hand=np.zeros(13,dtype=np.float64)
    pca=np.zeros(13,dtype=np.float64)
    ing=np.zeros(13,dtype=np.float64)
    
    grdSTATION.Initialized=False
    grdSTATION.saveIndex=(Ncohorts,Ntime,Nlarva,Nprey)
    grdSTATION.Ndays=Ntime/24.
    grdSTATION.Nlarva=Nlarva
    grdSTATION.Ncohorts=Ncohorts
    grdSTATION.Nprey=Nprey
    grdSTATION.larvaPsur=1.0
    grdSTATION.firstBatch=True
    
    t=1
    loopsDone = False
    grdSTATION.dayOfYear = -9
    print "RANDOM WEIGHT INITIALIZED"
    W[:,:,:,:]         = initWgt + ((initWgt/2.)*np.random.random_sample(W.shape))*randomWgt# milligram (5mm larvae) 
    W_AF[:,:,:,:]      = initWgt
    S[:,:,:,:]         = stomach_threshold*gut_size*W # 30% av max mageinnhold
    L[:,:,:,:]         = calculateLength(initWgt, 0.0)
    print W    
    larvaTime=[]
    released=[]
    for i in range(Ncohorts): released.append(0)
 
    depth=initDepth
    releasedCohorts=0
    isDead=[]
    for i in range(Ncohorts): isDead.append(False)
  
    """==================LOOP STARTS====================="""
    """Loop over all the days of interest"""
    
     
    for day in range(Ntime):
        if loopsDone : break
        """Save the julian date and the time index so that you can reset time between each chohort and indiviual you run"""
        oldJulian=julian
        oldT=t
        release=[];isReleased=[];releaseDate=[]
        
        if grdSTATION.seawifs is True:
            seawifsDate=grdSTATION.refDate + datetime.timedelta(seconds=julian)
            """Months goes from 1-12, so index needs to extract one"""
            index=int(seawifsDate.month) - 1
            maxSeawifsValue=seawifsArray[stationNumber,:].max()
            seawifsValue = seawifsArray[stationNumber,index]
            grdSTATION.relativeSeawifsValue = seawifsValue / maxSeawifsValue
            grdSTATION.seawifsValue = seawifsValue 
            #print "Food availability is defined by seawifs for month %i and station %s: %s"%(index+1,stationNumber,grdSTATION.relativeSeawifsValue)
           #TODO: add seawifs relative prey value to calculateGrowth function
            
        """loop over all the cohorts of interest"""
        noOneLeft=True
        for cohort in range(Ncohorts):
            """For each day, and cohort, loop over all the individuals"""
            a,b,c= checkReleased(grdSTATION.listOfReleaseDates[cohort],julian,grdSTATION)
            release.append(a); isReleased.append(b); releaseDate.append(c)
            
            if release[cohort] is True:
           #     print "\nReleasing a new cohort on %s"%(releaseDate[cohort])
                released[cohort]=t
                isReleased[cohort]=True
                releasedCohorts+=1
                noOneLeft=False
                
       # print "We have found %s cohorts that have been released"%(releasedCohorts)     
        for cohort in range(releasedCohorts):
            
            if isReleased[cohort] is False:
                continue
            elif released[cohort]!=0:
                activeCohorts, minNumberOfActiveCohort = findActiveCohorts(isDead,isReleased)
                
                for prey in range(Nprey):
                    
                    for ind in range(Nlarva):
                        if loopsDone : break
                        """We create an array that controls the entrance and exit points in time for
                        a given larvae"""
                        if startAndStopIndex[cohort,ind,0]<=missingValue:
                            startAndStopIndex[cohort,ind,0]=t-1
                            Age[cohort,ind,prey]=0
                      
                        isAliveBool, grdSTATION = isAlive(julian,Age[cohort,ind,prey],cohort,ind,t,prey,startAndStopIndex,isDead,grdSTATION)
                        if isAliveBool == False:
                            
                            if noOneLeft==True and ind==(Nlarva-1) and cohort==(releasedCohorts-1) and prey==(Nprey-1):
                                currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*24.)
                                delta=currentDate - grdSTATION.refDate 
                                julian=delta.days*86400 + delta.seconds
                                noOneLeft=True
                              #  print "Updating time as no cohort is alive currently",currentDate
                        else:
                            noOneLeft=False
                         
                            """We calculate swimspeed once every 24 hours and assume it does not change during that time-period"""
                            swimSpeed=0.261*(L[cohort,ind,t-1,prey]**(1.552*L[cohort,ind,t-1,prey]**(0.920-1.0)))-(5.289/L[cohort,ind,t-1,prey])
                            maxHourlyMove=((swimSpeed*(dt/1000.))/4.0)*deltaH # Divided by four compared to original IBM in fortran
                            maxHourlyMove = round(maxHourlyMove,lastDecimal)
                                
                            for hour in range(dt_per_day):
                                """Since we split hour into parts, the new hour variable is h"""   
                                h = deltaH*1.0*hour
                
                                oldDepth,optDepth,h_start,h_stop,NOptDepths,oldFitness = initBehavior(depth,maxHourlyMove)
                                gtime=grdSTATION.refDate + datetime.timedelta(seconds=julian)
                                
                                """LIGHT == Get the maximum light at the current latitude, and the current surface light at the current time.
                                These values are used in calculateGrowth and  predation.FishPredAndStarvation, but need
                                only be calcuated once per time step. """
                                tt = gtime.timetuple()
                                dayOfYear = float(tt[7]); month=float(tt[1]); hourOfDay = float(tt[3]); daysInYear=365.0
                                radfl0=0.0; maxLight=0.0; cawdir=0.0; clouds=0.0; sunHeight=0.0; surfaceLight=0.0
                                
                                radfl0,maxLight,cawdir = calclight.calclight.qsw(radfl0,maxLight,cawdir,clouds,grdSTATION.lat*np.pi/180.0,dayOfYear,daysInYear)
                                maxLight = maxLight/0.217 # Convert from W/m2 to umol/m2/s-1
                                sunHeight, surfaceLight = calclight.calclight.surlig(hourOfDay,maxLight,dayOfYear,grdSTATION.lat,sunHeight,surfaceLight)
                                """If you don't have fortran compiler:"""     
                                #surfaceLight = IOlight.surfaceLight(grdSTATION,julian,grdSTATION.lat,hour)
                                """Find the attenuation coefficient as a function of Chlorophyll-a values"""
                                if grdSTATION.seawifs is True:
                                    attCoeff=calculateAttenuationCoeff(grdSTATION.seawifsValue)
                                    beamAttCoeff=attCoeff*3.0
                                    #print "attenuation coeff: %s depth: %s chl: %s"%(attCoeff,depth,grdSTATION.relativeSeawifsValue)
                                        
                                """OPTIMAL DEPTH == Find the optimal depth to stay at given ratio of ingestion and mortality rate"""
                                for findDepth in range(NOptDepths+1):
                                    """Current depth:"""
                                    depth=h_start+findDepth*deltaZ
                                    """Light level at current depth:"""
                                    Eb = surfaceLight*np.exp(attCoeff*(-depth))
                                    """Calculate growth and ingestion at the current depth:"""
                                    
                                    ing, GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness,Ndata = calculateGrowth(julian,Eb,deltaH,depth,h,grdSTATION,dt,W[cohort,ind,t-1,prey],
                                                                                                         L[cohort,ind,t-1,prey],julianIndex,
                                                                                                         julianFileA,julianFileB,R,enc,hand,ing,pca,
                                                                                                         S[cohort,ind,t-1,prey],prey)
                         
                                    """Calculate mortality at the current depth layer"""
                                    mortality, didStarve, dead = predation.FishPredAndStarvation(grdSTATION,deltaH,FishDens,L[cohort,ind,t-1,prey]*mm2m,W[cohort,ind,t-1,prey],
                                                                                attCoeff,Eb,dt,ing,stomachFullness)
                                  
                                    """Calculate fitness at the current depth layer"""
                                    F = sum(ing)/W[cohort,ind,t-1,prey]

                                    
                                    optDepth,oldFitness = getBehavior(stomachFullness,
                                                                 F,mortality,L[cohort,ind,t-1,prey],
                                                                 oldFitness,depth,optDepth)
                                   
                                """Set the optimal depth equal to result of optimal loop and update light at depth"""    
                                depth=optDepth
                                if grdSTATION.seawifs is True:
                                    attCoeff=calculateAttenuationCoeff(grdSTATION.seawifsValue)
                                    beamAttCoeff=attCoeff*3.0
                                        
                                Eb = surfaceLight*np.exp(attCoeff*(-depth))
                                """Now recalculate the growth and mortality at the optimal depth layer"""
                                ing, GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness,Ndata = calculateGrowth(julian,Eb,deltaH,depth,h,grdSTATION,dt,W[cohort,ind,t-1,prey],
                                                                                                     L[cohort,ind,t-1,prey],julianIndex,
                                                                                                     julianFileA,julianFileB,R,enc,hand,ing,pca,
                                                                                                     S[cohort,ind,t-1,prey],prey)
                                
                                mortality, didStarve, dead = predation.FishPredAndStarvation(grdSTATION,deltaH,FishDens,L[cohort,ind,t-1,prey]*mm2m,
                                                                                       W[cohort,ind,t-1,prey],attCoeff,
                                                                                       Eb,dt,ing,stomachFullness)
                               
                                S[cohort,ind,t,prey] = max(0.0, min(gut_size*W[cohort,ind,t-1,prey],S[cohort,ind,t-1,prey] + sum(ing[:])))
                                ingrate[cohort,ind,t,prey] = max(0.0, (S[cohort,ind,t,prey] - S[cohort,ind,t-1,prey])/W[cohort,ind,t-1,prey])
                              
                                W[cohort,ind,t,prey] = W[cohort,ind,t-1,prey] + min(GR_gram + meta,S[cohort,ind,t,prey]*assi) - meta 
                                S[cohort,ind,t,prey] = max(0.0, S[cohort,ind,t,prey] - ((W[cohort,ind,t,prey] - W[cohort,ind,t-1,prey]) + meta)/assi) # + 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta)/assi;
                                L[cohort,ind,t,prey] = calculateLength(W[cohort,ind,t,prey], L[cohort,ind,t-1,prey])
                               
                                W_AF[cohort,ind,t,prey] = W_AF[cohort,ind,t-1,prey] + (np.exp(GR)-1)*W_AF[cohort,ind,t-1,prey]
                                larvaTdata[cohort,ind,t,prey] = Tdata
                                larvaDepth[cohort,ind,t,prey] = depth
                                larvaNauplii[cohort,ind,t,prey] = Ndata
                                
                                if didStarve > 1: larvaPsur[cohort,ind,t,prey]=0.0
                                else:larvaPsur[cohort,ind,t,prey]=larvaPsur[cohort,ind,t-1,prey]*(np.exp(-mortality))
                                #
                                #print 'Prob.survival at depth : %3.2f for time %s => %3.6f for size %s'%(depth,
                                #                                                            grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                                #                                                            ,larvaPsur[cohort,ind,t,prey]*100.,
                                #                                                            L[cohort,ind,t,prey])
                                                                                            
                                SGR[cohort,ind,t,prey]=((W[cohort,ind,t,prey]-W[cohort,ind,t-1,prey])/W[cohort,ind,t-1,prey])*100.0
                                
                                if ind==0 and cohort==minNumberOfActiveCohort and prey==0: #releasedCohorts-1:
                                    larvaTime.append(julian/3600.)
                                    
                                julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION,clim)
                                   
                                s1=int(startAndStopIndex[cohort,ind,0])
                                
                                Age[cohort,ind,prey] =julian/3600.-larvaTime[s1]
                                
                                """Update time"""
                                currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                                delta=currentDate - grdSTATION.refDate
                                 
                                julian=delta.days*86400 + delta.seconds
                                t +=1
                             
                            now=grdSTATION.refDate + datetime.timedelta(seconds=julian)
                            last=grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive)
                  
                            if now == grdSTATION.endDate or now == last:
                               
                                if ind==(Nlarva-1) and cohort==(releasedCohorts-1) and prey==(Nprey-1):
                                    message='---> Finished IBMtime on date: %s. Writing results to file....'%(grdSTATION.refDate + datetime.timedelta(seconds=julian))
                                    showProgress(Ntime,Ntime,message)
                                    larvaTime.append(julian/3600.)

                                    startAndStopIndex[cohort,:,1]=t
                                    IOwrite.writeStationFile(deltaH,deltaZ,grdSTATION,larvaTime,W,L,SGR,larvaTdata,larvaNauplii,larvaDepth,W_AF,larvaPsur,outputFile,startAndStopIndex)
                                    move(outputFile, "results/"+outputFile)
                                    loopsDone = True
                                    
                                    print "Result file is saved as %s"%("results/"+outputFile)
                                    break
                                    
                            if ind==(Nlarva-1) and cohort==(releasedCohorts-1) and prey==(Nprey-1):
                                julian=julian
                                t=t
                            else:
                                julian=oldJulian
                                t=oldT
                            julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION,clim)
                            hour=0
                    ind=0
        if loopsDone is False:
            """Show progress indicator"""
            message='---> running IBMtime-step %s of %s with %s released cohorts (%s finished)'%(t,grdSTATION.refDate + datetime.timedelta(seconds=julian),releasedCohorts,minNumberOfActiveCohort)
            showProgress(t,Ntime,message)
    
if __name__=="__main__":
   
    try:
        import psyco
        psyco.log()
        psyco.profile()
    except ImportError:
        pass
    
    stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_NorthSea.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_Iceland.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_Lofoten.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc"]
    """Notice that if you use seawifs data, the order of the stations here have to match
    the order of station seawifs data in seawifs file (see extractSeaWifs.py file)"""
    stationNames=['North Sea', 'Iceland','Lofoten', 'Georges Bank']
    
    stationNumber=0
    
    events = ['COLD', 'WARM']
    #event = 'REGULAR RUN'
    
    if events[0] !='COLD':
        for station, stationName in zip(stations,stationNames):
            
            ibm(station, stationName, stationNumber, event)
            stationNumber+=1
    else:
        for station, stationName in zip(stations,stationNames):
            for event in events:
                ibm(station, stationName, stationNumber, event)
            
            stationNumber+=1
