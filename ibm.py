import os, sys, string
from shutil import move
import netCDF4
from netCDF4 import num2date, Dataset
import types
import numpy as np
import IOwrite
import datetime as datetime
from pylab import *
import help

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
import IOnetcdf
import predation
from initLarva import *
"""Import f2py Fortran modules"""
import calclight
import perception
import bioenergetics
import IOstation
from memoize import memoize

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
    Finish      = False
    useAverageFile = False
    fileNameIn  = station

    """Option for using seawifs data as prey abundance"""
    seawifs=False
    seawifsFileName="/Users/trond/Projects/seawifs/chlo-stations.nc"
    if seawifs==True:
        seawifsArray=IOnetcdf.getSeaWifs(seawifsFileName)
    else:
        seawifsArray=None

    """Open the netcdf file if it existst."""
    cdf = Dataset(fileNameIn)

    """Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function."""
    grdSTATION = grd.grdClass(fileNameIn,"STATION")

    if useAverageFile==True:
        inFile="/Users/trond/Projects/arcwarm/SODA/soda2average/clim/averageSODA1961-1990.nc"
        aveTemp = IOnetcdf.getAverageValues(inFile,grdSTATION.lon,grdSTATION.lat)
        grdSTATION.aveT=aveTemp

    if event == "REGULAR RUN":
        varlist=['temp','salt','u','v','taux','tauy'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
        startDate = datetime.datetime(1970,4,1,0,0,0)
        endDate   = datetime.datetime(1970,4,6,0,0,0)

    if event == "CLIM RUN":
        varlist=['temp','salt','u','v'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
        startDate = datetime.datetime(1948,1,15,1,0,0)
        endDate   = datetime.datetime(1948,12,14,0,0,0)

    if event == "ESM RUN":
        varlist=['temp','salt','nlg','nsm','chla','taux','tauy'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
        startDate = datetime.datetime(2002,4,15,1,0,0)
        endDate   = datetime.datetime(2002,8,14,0,0,0)
        grdSTATION.chlaValue=0.0

    if event == 'COLD' or event=='WARM':
        varlist=['temp','salt','u','v','taux','tauy']
        startDate, endDate = coldWarmEvent(stationName,event)

    """Get the time information and find the indices for start and stop data to extract relative to
    the time period wanted. Takes input data:
    startDate="DD/MM/YYYY" and endDate="DD/MM/YYYY" or if none given, finds all date"""
    getTimeIndices(cdf,grdSTATION,startDate,endDate,event)

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
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,log,stationName)
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

    for i in xrange(NReleaseDatesInit):
        date=startDate+datetime.timedelta(daysBetweenReleases*(i+1))
        if endDate >= date + datetime.timedelta(days=NDaysAlive):
            listOfReleaseDates.append(date)

    NReleaseDates=len(listOfReleaseDates)
    grdSTATION.listOfReleaseDates=listOfReleaseDates

    print "\nThis simulation will release a total of %s cohorts on the following dates:"%(NReleaseDates)
    for h in xrange(NReleaseDates):
        print "--> %s"%(listOfReleaseDates[h])
    print "\n"

    print "IBM will run from date %s to %s\n"%(grdSTATION.listOfReleaseDates[0] + datetime.timedelta(days=0),
                                              grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive))
    print "-----------------------------------------------------------------"

    grdSTATION.seawifs=seawifs
    grdSTATION.stationNumber=0

    infoOnResolution(grdSTATION)

    return grdSTATION, outputFile, NReleaseDates, Finish, seawifsArray

def findActiveCohorts(isDead,isReleased):
    co=0; co2=0
    for check in xrange(len(isDead)):
        if isDead[check]==True: co+=1

    for check in xrange(len(isReleased)):
        if isReleased[check]==True: co2+=1
    activeCohorts=np.abs(co2-co)
    minNumberOfActiveCohort = co
    #print "209: There are currently %s active cohorts in the simulation with min cohort %s"%(activeCohorts,minNumberOfActiveCohort)

    return minNumberOfActiveCohort

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

    #print "236:Will release:",release, " Has been released",isReleased, " at date:",releaseDate," current date:", grdSTATION.refDate + datetime.timedelta(seconds=julian)
    #print "237:Current date:", grdSTATION.refDate + datetime.timedelta(seconds=julian)
    return release, isReleased, releaseDate


def isAlive(julian,larvaAge,cohort,isDead):
    if isDead[cohort]==False:
        #print "Larva is %s days old and belongs to cohort %s %s"%(larvaAge/(24.),cohort, NDaysAlive)

        if larvaAge/(24.) >= int(NDaysAlive):
            isAliveBool=False
            isDead[cohort]=True
        else:
            isAliveBool=True
    else:
        isAliveBool=False

    return isAliveBool

@memoize
def getDepthIndex(grdSTATION,depth):
    depthFOUND=False; d=0
    
    while d < len(grdSTATION.depth) and depthFOUND==False:

        if abs(grdSTATION.depth[d]) <= abs(depth) <= abs(grdSTATION.depth[d+1]) and depthFOUND==False:
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

def initBehavior(depth,maxHourlyMove,grdSTATION):
    oldDepth=depth  #The previous depth for a given larvae for a time step inside optimal loop
    optDepth=depth  # The optimal depth for a given larvae for a time step
    h_start, h_stop = getMaxMinDepths(oldDepth,maxHourlyMove,grdSTATION)
    NOptDepths=int((abs(oldDepth-h_start) + abs(h_stop-oldDepth) )/float(deltaZ))

    sampleDepths=[]
    maxDiff=0.0
    if oldDepth  == 0.0:
        for i in xrange(NOptDepths+1):
            newDepth = oldDepth + float(deltaZ)*(float(i))
            if newDepth <0: newDepth=0.0
            if newDepth < h_start: newDepth=h_start
            if newDepth > h_stop: newDepth=h_stop
            if newDepth not in sampleDepths:
                sampleDepths.append(newDepth)
                if abs(newDepth-oldDepth) > maxDiff: maxDiff=abs(newDepth-oldDepth)
    else:
        # Move down into the water column limited by h_stop
        for i in xrange(NOptDepths/2):
            newDepth = oldDepth + float(deltaZ)*(NOptDepths/2 - float(i))
            if newDepth <0: newDepth=0.0

            if newDepth > h_stop: newDepth=h_stop
            if newDepth not in sampleDepths:
                sampleDepths.append(newDepth)
                if abs(newDepth-oldDepth) > maxDiff: maxDiff=abs(newDepth-oldDepth)
        for i in xrange(NOptDepths/2+1):
            newDepth = oldDepth - float(deltaZ)*((float(i)))
            if newDepth <0: newDepth=0.0
            if newDepth < h_start: newDepth=h_start
            if newDepth not in sampleDepths:
                sampleDepths.append(newDepth)
                if abs(newDepth-oldDepth) > maxDiff: maxDiff=abs(newDepth-oldDepth)

    #for j in range(len(sampleDepths)):
    #  print "depth %3.3f diffZ %3.3f"%(sampleDepths[j], abs(oldDepth-sampleDepths[j]))
    fitness=-9999.9

    return oldDepth,optDepth,sampleDepths,maxDiff,fitness

def getBehavior(stomachFullness,F,m,length,oldFitness,depth,optDepth):
    """Rule 4 of Behavioral Ecology paper - Kristiansen et al. 2009"""
    T=min(0.9,0.3+1000.0*(1+(length)*np.exp(length))**(-1))

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


def calculateLength(Larval_wgt, Larval_m):
    oldLarval_m = Larval_m
    Larval_m = np.exp(2.296 + 0.277*np.log(Larval_wgt) - 0.005128*np.log(Larval_wgt)**2)
    if Larval_m < oldLarval_m:
        Larval_m = oldLarval_m
        #print 'Kept the old length', oldLarval_m, Larval_m
    return Larval_m


def calculateGrowth(Tdata,NLGdata,windX,windY,julian,Eb,deltaH,depth,hour,grdSTATION,dt,Larval_wgt, Larval_mm,
                    R,enc,hand,ing,pca,Spre,prey,event):

    """Mouthsize of larvae. The larvae can only capture calanus of sizes less
    than the larval mouthsize. Folkvord et al."""
    m = np.exp (-3.27+1.818*np.log(Larval_mm)-0.1219*(np.log(Larval_mm))**2.)
    """Calculate metabolism"""
    meta = dt*2.38e-7*np.exp(0.088*Tdata)*((Larval_wgt*mg2ug)**(0.9)*0.001)*deltaH
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
    #currentDate = grdSTATION.refDate + datetime.timedelta(seconds=julian)
    #if event!="ESM RUN":
    #    Tanomaly = (Tdata - grdSTATION.aveT[int(currentDate.month-1)]) / (grdSTATION.aveT.max() - grdSTATION.aveT.min())
    #else:
    Tanomaly = 0.0

    #if grdSTATION.seawifs is True:
    #    if grdSTATION.relativeSeawifsValue > 0.01:
    #        P = Tanomaly + grdSTATION.relativeSeawifsValue
    #    else:
    #        P = grdSTATION.relativeSeawifsValue
    #    P = max(0.0,P)
    if event=="ESM RUN":
        """Use large phytoplankton as proxy for zooplankton"""
        P=0.45 * NLGdata * 1.e4

    """VISUAL == PERCEPTION of PREY calculations==============================="""
    Em = (Larval_mm**2.0)/(contrast*0.1*0.2*0.75) #Size-specific sensitivity of the visual system (Fiksen,02)
    ing[:]=0.0; enc[:]=0.0; hand[:]=0.0; pca[:]=0.0

    for j in xrange(16):

        IER=0; R[j]=0.0 #R[j]=np.sqrt(Em*contrast*(prey_AREA[j])*(Eb/(Ke_larvae+Eb)))
        """All input to getr is either in m (or per m), or in mm (or per mm)"""
        R[j], IER = perception.perception.getr(R[j],beamAttCoeff/m2mm,contrast,prey_AREA[j],Em,Ke_larvae,Eb, IER)
        pca[j] = max(0.0,min(1.0,-16.7*(prey_LENGTH[j]/Larval_mm) + 3.0/2.0))
        """If you don't have fortran compiler use python routine for visual range"""
        #visual = np.sqrt(Em*contrast*(Ap_calanus*m2mm)*(Eb/(Ke_larvae+Eb)))
        #R[j] = (IOlight.getPerceptionDistance(Em,attCoeff,Ke_larvae,Ap_calanus,Eb))*m2mm

    """Calculate turbulence based on wind stress"""
    epsilon=(5.82*1.E-9*((np.sqrt(windX**2 + windY**2)))**3.)/(depth+0.1)
    omega = 1.9*(epsilon*R[:]*mm2m)**0.667
    omega =  omega * m2mm # From m/s to mm/s

    """Calculate handling time, encounter rate, and probability of capture"""
    hand[:] = 0.264*10**(7.0151*(prey_LENGTH[:]/Larval_mm)) # Walton 1992
    enc[:] = ((0.667*np.pi*(R[:]**3.)*f + np.pi*(R[:]**2.)*np.sqrt(prey_LENGTH[:]**2.+ 2.*omega**2.)*f*tau) * (prey_D[:]*((prey+1)*MultiplyPrey*P))* ltr2mm3)

    """Calculate ingestion rate"""
    ing[:] = (dt*enc[:]*pca[:]*prey_WGT[:]*micro2m / (1 + hand[:]))*deltaH

    for j in xrange(16):
        print "j:",j, " enc:",enc[j]*dt*deltaH, " mm:",Larval_mm, "pca:",pca[j], "SPpl:",prey_D[j]
        print "j:",j, " ing:",ing[j], "visual:",R[j]
    """Calculate stomach fullness"""
    stomachFullness =  (min(gut_size*Larval_wgt,Spre + sum(ing)))/(Larval_wgt*gut_size)
    
    return ing,GR_gram,GR,meta,assi,stomachFullness

def getMaxMinDepths(oldDepth,maxHourlyMove,grdSTATION):
    """The value you give deepstDepthAllowed will be the absolute deepest
    depth the larvae will move to. Useful to set this depth equal to bottom.
    grdSTATION.deepestDepthAllowed is defined from data in function IOnetcdf.getData"""

    if (oldDepth>maxHourlyMove and oldDepth+maxHourlyMove < grdSTATION.deepestDepthAllowed):
        h_start = oldDepth - maxHourlyMove
        h_end   = oldDepth + maxHourlyMove

    elif (oldDepth <= maxHourlyMove):
        h_start=0.0
        h_end=min(grdSTATION.deepestDepthAllowed,oldDepth + maxHourlyMove)
    elif (oldDepth + maxHourlyMove >= grdSTATION.deepestDepthAllowed):
        h_start=max(1,oldDepth-maxHourlyMove)
        h_end=grdSTATION.deepestDepthAllowed
    else:
        print 'STOP: Depth error (ibm.py:getMaxMinDepths)'

    return h_start, h_end


def calculateAttenuationCoeff(Chla):
    """See Fiksen et al. 2002 for details (relationship taken from Riley 1956)
    units of Chla is mg/m3"""
    k0 = 0.1
    return k0 + 0.054*Chla**(2./3.) + 0.0088*Chla

def getTimeIndices(cdf,grdSTATION,startDate=None,endDate=None,event=None):
    """Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    if event=="ESM RUN":
        refDate    = datetime.datetime(2001,1,1,0,0,0)
    else:
        refDate    = datetime.datetime(1948,1,1,0,0,0)
    grdSTATION.refDate = refDate

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

        for i in xrange(grdSTATION.time.shape[0]):

            if grdSTATION.time[i]*86400  < JDstart:
                FOUND=False
                continue
            elif grdSTATION.time[i]*86400  >= JDstart and FOUND is False:
                print "\nFound first time index that fits start point:"
                print "%s at index %i => %s"%(refDate + datetime.timedelta(seconds=grdSTATION.time[i-1]*86400),i-1,JDstart)
                grdSTATION.startIndex=i-1
                grdSTATION.start_date = refDate + datetime.timedelta(seconds=grdSTATION.time[i-1]*86400)
                FOUND=True

        for i in xrange(grdSTATION.time.shape[0]):
            if grdSTATION.time[i]*86400 < JDend:
                FOUND=False
                continue
            elif grdSTATION.time[i]*86400 >= JDend and FOUND is False:
                print "Found first time index that fits end point:"
                print "%s at index %i => %s\n"%(refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400),i,JDend)
                grdSTATION.endIndex=i

                grdSTATION.end_date = refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400)
                FOUND=True

    print "Total number of days extracted: %s"%(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])

    grdSTATION.JDend=JDend
    grdSTATION.JDstart=JDstart
    grdSTATION.totalDays=(grdSTATION.time[grdSTATION.endIndex]-grdSTATION.time[grdSTATION.startIndex])


def updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION):
    """This routine makes sure that we are reading at the correct time indices
    in the input data. If we have moved forward in time enough to move past the current time window of two
    time stamps, then update the indices and give the new indices.
    Remember that input time format from netCDF files must be
    'Days since 1948.1.1'
    """
    refDate    = datetime.datetime(1948,1,1,0,0,0)
   # print "searching: \n%s \n%s \n%s"%(refDate + datetime.timedelta(seconds=julianFileA),
   #                                    refDate + datetime.timedelta(seconds=julian),
   #                                    refDate + datetime.timedelta(seconds=julianFileB))

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
    julianFileA=grdSTATION.time[julianIndex]*86400
    julianFileB=grdSTATION.time[julianIndex+1]*86400

    return julianFileA, julianFileB, julianIndex, Finish


def showProgress(t,Ntime,message):
    p=int( ((t*1.0)/(1.0*(Ntime-24*dt_per_day)))*100.)
    progress.render(p,message)

def ibm(station,stationName,stationNumber,event):

    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION, outputFile, Ncohorts, Finish, seawifsArray = init(station,stationName,event)

    """Start time in days, hours, and seconds since 1948,1,1"""
    time=grdSTATION.startDate - grdSTATION.refDate

    """ julian is the time used to control the larvae, while
    julianFileA and julianFileB are the time stamps of inbetween
    the larvae recides relative to time from file. Temperature is
    interpolated to date julian from julianFileA and julianFileB"""
    julian=time.days*86400. + time.seconds

    julianFileA=grdSTATION.time[0]*86400.
    julianFileB=grdSTATION.time[1]*86400.
    Ntime=int(grdSTATION.delta) # Number of time-steps to use in arrays

    julianIndex=0

    """Create arrays to store all information in"""
    Tlarva=(NDaysAlive+1)*(1.0*dt_per_day) + 1

    print "Each larva will store information for %s timesteps: deltaH=%s"%(Tlarva,Tlarva/dt_per_day)
    W            =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    larvaPsur    =(np.ones((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))
    larvaDepth   =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    W_AF         =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    S            =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    Age          =(np.ones((Ncohorts,Nlarva,Nprey),dtype=np.float64))*missingValue
    startAndStop =(np.ones((Ncohorts,Nlarva,2),dtype=np.float64))*missingValue
    L            =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue

    ingrate      =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    Psurvive     =(np.ones((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    SGR          =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    larvaTdata   =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue
    larvaNauplii =(np.zeros((Ncohorts,Nlarva,Tlarva,Nprey),dtype=np.float64))*missingValue

    R=np.zeros(16,dtype=np.float64)
    enc=np.zeros(16,dtype=np.float64)
    hand=np.zeros(16,dtype=np.float64)
    pca=np.zeros(16,dtype=np.float64)
    ing=np.zeros(16,dtype=np.float64)

    grdSTATION.Initialized=False
    #grdSTATION.saveIndex=(Ncohorts,Ntime,Nlarva,Nprey)
    grdSTATION.Ndays=Ntime/24.
    grdSTATION.Nlarva=Nlarva
    grdSTATION.Ncohorts=Ncohorts
    grdSTATION.Nprey=Nprey
    grdSTATION.larvaPsur=1.0
    grdSTATION.firstBatch=True
    grdSTATION.firstRun=True

    t=1
    loopsDone = False
    grdSTATION.dayOfYear = -9
    if randomWgt==1: print "RANDOM WEIGHT INITIALIZED"

    W[:,:,:,:]         = initWgt + ((initWgt*0.2)* (np.random.random_sample(W.shape))*np.random.random_integers(-1,1,W.shape)*randomWgt) # milligram (5mm larvae)
    W_AF[:,:,:,:]      = initWgt
    S[:,:,:,:]         = stomach_threshold*gut_size*W # 30% av max mageinnhold
    L[:,:,:,:]         = calculateLength(initWgt, 0.0)

    larvaTime=[]

    depth=initDepth
    releasedCohorts=0
    isDead=[]
    minNumberOfActiveCohort=0
    releasedCohorts=0
    released=[]
    for i in xrange(Ncohorts):  released.append(0)
   
    for i in xrange(Ncohorts): isDead.append(False)

    """==================LOOP STARTS====================="""
    """Loop over all the days of interest"""

    for day in xrange(Ntime):
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

        """loop over all the cohorts of interest"""
        noOneLeft=True
        currentCohort=0
        
        for cohort in xrange(minNumberOfActiveCohort,releasedCohorts+1):
           
            """For each day, and cohort, loop over all the individuals"""
            if cohort >= len(grdSTATION.listOfReleaseDates):
                continue
            else:
                a,b,c= checkReleased(grdSTATION.listOfReleaseDates[cohort],julian,grdSTATION)
                release.append(a); isReleased.append(b); releaseDate.append(c); 
                
                if release[currentCohort] is True:
               #     print "\nReleasing a new cohort on %s"%(releaseDate[currentCohort])
                    released[cohort]=t
                    isReleased[currentCohort]=True
                   
                    releasedCohorts+=1
                    noOneLeft=False
                
                currentCohort+=1
            
   #     print "We have found %s cohorts that have been released"%(releasedCohorts)
       
        currentCohort=0
        for cohort in xrange(minNumberOfActiveCohort,releasedCohorts):
            #print "Checking # %s : cohorts from %s to %s : currentCohort %s"%(cohort,minNumberOfActiveCohort,releasedCohorts,currentCohort)
          
            if isReleased[currentCohort] is False:
                continue
            elif released[cohort]!=0:
                minNumberOfActiveCohort = findActiveCohorts(isDead,isReleased)
        
                for prey in xrange(Nprey):

                    for ind in xrange(Nlarva):
                        if day==0:
                            depth =initDepth * np.random.random_sample()
                            print "initdepth is %s and depth is %s"%(initDepth,depth)
                        if loopsDone : break
                        """We create an array that controls the entrance and exit points in time for
                        a given larvae"""

                        if startAndStop[cohort,ind,0]<=missingValue:
                            startAndStop[cohort,ind,0]=t-1
                            Age[cohort,ind,prey]=0
                        """This is the index that counts the time for each individual larva"""
                        larvaIndex=t-startAndStop[cohort,ind,0]

                        isAliveBool = isAlive(julian,Age[cohort,ind,prey],cohort,isDead)
                        if isAliveBool == False:

                            if noOneLeft==True and ind==(Nlarva-1) and cohort==(releasedCohorts-1) and prey==(Nprey-1):
                                currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*24.)
                                delta=currentDate - grdSTATION.refDate
                                julian=delta.days*86400 + delta.seconds
                                noOneLeft=True
                           #     print "Updating time as no cohort is alive currently",currentDate
                        else:
                            noOneLeft=False
                            #print "737: Larval index t=%s, t-1=%s and t=%s"%(larvaIndex, larvaIndex-1, t)
                            """We calculate swimspeed once every 24 hours and assume it does not change during that time-period"""
                            swimSpeed=0.261*(L[cohort,ind,larvaIndex-1,prey]**(1.552*L[cohort,ind,larvaIndex-1,prey]**(0.920-1.0)))-(5.289/L[cohort,ind,larvaIndex-1,prey])
                            maxHourlyMove=((swimSpeed*(dt/1000.))/2.0)*deltaH # Divided by four compared to original IBM in fortran
                            maxHourlyMove = round(maxHourlyMove,lastDecimal)

                            for hour in xrange(dt_per_day):

                                """Since we split hour into parts, the new hour variable is h"""
                                h = deltaH*1.0*hour

                                oldDepth,optDepth,sampleDepths,maxDiffZ,oldFitness = initBehavior(depth,maxHourlyMove,grdSTATION)
                                gtime=grdSTATION.refDate + datetime.timedelta(seconds=julian)

                                """LIGHT == Get the maximum light at the current latitude, and the current surface light at the current time.
                                These values are used in calculateGrowth and  predation.FishPredAndStarvation, but need
                                only be calcuated once per time step. """
                                tt = gtime.timetuple()
                                dayOfYear = float(tt[7]); month=float(tt[1]); hourOfDay = float(tt[3]); daysInYear=365.0
                                radfl0=0.0; maxLight=0.0; cawdir=0.0; clouds=0.0; sunHeight=0.0; surfaceLight=0.0

                                """Calculate the maximum and average light values for a given geographic position
                                for a given time of year. Notcie we use radfl0 instead of maximum values maxLight
                                in light caclualtions. Seemed better to use average values than extreme values.
                                NOTE: Convert from W/m2 to umol/m2/s-1"""
                                radfl0,maxLight,cawdir = calclight.calclight.qsw(radfl0,maxLight,cawdir,clouds,grdSTATION.lat*np.pi/180.0,dayOfYear,daysInYear)
                                maxLight = radfl0/0.217
                                sunHeight, surfaceLight = calclight.calclight.surlig(hourOfDay,maxLight,dayOfYear,grdSTATION.lat,sunHeight,surfaceLight)

                                """Find the attenuation coefficient as a function of Chlorophyll-a values"""
                                if grdSTATION.seawifs is True:
                                    attCoeff=calculateAttenuationCoeff(grdSTATION.seawifsValue)
                                    beamAttCoeff=attCoeff*3.0
                                    #print "attenuation coeff: %s depth: %s chl: %s"%(attCoeff,depth,grdSTATION.relativeSeawifsValue)
                                if event=="ESM RUN":
                                    attCoeff=calculateAttenuationCoeff(grdSTATION.chlaValue)
                                   # print grdSTATION.chlaValue, attCoeff
                                    beamAttCoeff=attCoeff*3.0

                                """OPTIMAL DEPTH == Find the optimal depth to stay at given ratio of ingestion and mortality rate"""
                                for depth in sampleDepths:

                                    diffZ=abs(oldDepth - depth)

                                    """Light level at current depth:"""
                                    Eb = surfaceLight*np.exp(attCoeff*(-depth))
                                    """Calculate growth and ingestion at the current depth:"""
                                    
                                    """==>INDEX calculations==============================="""
                                    depthIndex1, depthIndex2, dz1, dz2 = getDepthIndex(grdSTATION,depth)
    
                                    if event != "ESM RUN":
                                        Tdata,Sdata,Udata,Vdata,windX,windY = IOnetcdf.getData(julian,julianIndex,julianFileA,
                                                                          julianFileB,dz1,dz2,depthIndex1,
                                                                          depthIndex2,grdSTATION,event)
                                    else:
                                        Tdata,NSMdata,NLGdata,CHLAdata,windX,windY = IOnetcdf.getData(julian,julianIndex,julianFileA,
                                                                          julianFileB,dz1,dz2,depthIndex1,
                                                                          depthIndex2,grdSTATION,event)
                                        grdSTATION.chlaValue=CHLAdata

                                    print "Calling calulateGrowth"
                                    ing, GR_gram,GR,meta,assi,stomachFullness= calculateGrowth(Tdata,NLGdata,windX,windY,julian,Eb,deltaH,depth,h,grdSTATION,dt,
                                                                                                W[cohort,ind,larvaIndex-1,prey],
                                                                                                L[cohort,ind,larvaIndex-1,prey],R,enc,hand,ing,pca,
                                                                                                S[cohort,ind,larvaIndex-1,prey],prey,event)
                                    print "Calling bioenergetics.growth"
                                    II=len(prey_LENGTH)
                                    print Tdata,NLGdata,windX,windY
                                    print S[cohort,ind,larvaIndex-1,prey],prey_AREA,prey_LENGTH
                                    print prey_D,prey_WGT,L[cohort,ind,larvaIndex-1,prey],W[cohort,ind,larvaIndex-1,prey]
                                    print sec2day, mg2ug,ltr2mm3, m2mm,mm2m
                                    print f,tau,micro2m,Eb,contrast,Ke_larvae
                                    print beamAttCoeff,deltaH,dt,II,depth,prey,gut_size
                                    print "size",II
                           
                                    ing,GR_gram,GR,meta,assi,stomachFullness,zoop = bioenergetics.bioenergetics.growth(Tdata,NLGdata,windX,windY,
                                                                                    S[cohort,ind,larvaIndex-1,prey],prey_AREA,prey_LENGTH,
                                                                                    prey_D,prey_WGT,L[cohort,ind,larvaIndex-1,prey],W[cohort,ind,larvaIndex-1,prey],
                                                                                    sec2day, mg2ug,ltr2mm3, m2mm,mm2m,f,tau,micro2m,Eb,contrast,Ke_larvae,
                                                                                    beamAttCoeff,deltaH,dt,II,depth,prey,gut_size)
                                    

                                    """Calculate mortality at the current depth layer"""
                                    mortality, didStarve, dead = predation.FishPredAndStarvation(grdSTATION,deltaH,FishDens,
                                                                                                 L[cohort,ind,larvaIndex-1,prey]*mm2m,W[cohort,ind,larvaIndex-1,prey],
                                                                                                 attCoeff,Eb,dt,ing,stomachFullness)

                                    """Calculate the cost of being active. The more you swim the more you use energy in raltive ratio to
                                    routine metabolims. Trond Kristiansen, 23.03.2010"""
                                    activityCost= (diffZ/maxDiffZ)*(meta*costRateOfMetabolism)


                                    """Calculate fitness at the current depth layer"""
                                    F = (sum(ing)-activityCost)/W[cohort,ind,larvaIndex-1,prey]


                                    optDepth,oldFitness = getBehavior(stomachFullness,
                                                                 F,mortality,L[cohort,ind,larvaIndex-1,prey],
                                                                 oldFitness,depth,optDepth)

                                """Set the optimal depth equal to result of optimal loop and update light at depth"""
                                diffZ=abs(depth-optDepth)
                                depth=optDepth

                                Eb = surfaceLight*np.exp(attCoeff*(-depth))
                                """Now recalculate the growth and mortality at the optimal depth layer"""
                                ing, GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness,Ndata = calculateGrowth(julian,Eb,deltaH,depth,h,grdSTATION,dt,W[cohort,ind,larvaIndex-1,prey],
                                                                                                     L[cohort,ind,larvaIndex-1,prey],julianIndex,
                                                                                                     julianFileA,julianFileB,R,enc,hand,ing,pca,
                                                                                                     S[cohort,ind,larvaIndex-1,prey],prey,event)

                                mortality, didStarve, dead = predation.FishPredAndStarvation(grdSTATION,deltaH,FishDens,L[cohort,ind,larvaIndex-1,prey]*mm2m,
                                                                                       W[cohort,ind,larvaIndex-1,prey],attCoeff,
                                                                                       Eb,dt,ing,stomachFullness)

                                """Calculate the cost of the activity performed during last hour for vertical behavior"""
                                activityCost= (diffZ/maxDiffZ)*(meta*costRateOfMetabolism)

                                S[cohort,ind,larvaIndex,prey] = max(0.0, min(gut_size*W[cohort,ind,larvaIndex-1,prey],S[cohort,ind,larvaIndex-1,prey] + sum(ing[:])))
                                ingrate[cohort,ind,larvaIndex,prey] = max(0.0, (S[cohort,ind,larvaIndex,prey] - S[cohort,ind,larvaIndex-1,prey])/W[cohort,ind,larvaIndex-1,prey])

                                W[cohort,ind,larvaIndex,prey] = W[cohort,ind,larvaIndex-1,prey] + min(GR_gram + meta,S[cohort,ind,larvaIndex,prey]*assi) - meta - activityCost

                                S[cohort,ind,larvaIndex,prey] = max(0.0, S[cohort,ind,larvaIndex,prey] - ((W[cohort,ind,larvaIndex,prey] - W[cohort,ind,larvaIndex-1,prey]) + meta)/assi)
                                L[cohort,ind,larvaIndex,prey] = calculateLength(W[cohort,ind,larvaIndex,prey], L[cohort,ind,larvaIndex-1,prey])

                                W_AF[cohort,ind,larvaIndex,prey] = W_AF[cohort,ind,larvaIndex-1,prey] + (np.exp(GR)-1)*W_AF[cohort,ind,larvaIndex-1,prey]
                                larvaTdata[cohort,ind,larvaIndex,prey] = Tdata
                                larvaDepth[cohort,ind,larvaIndex,prey] = depth
                                larvaNauplii[cohort,ind,larvaIndex,prey] = Ndata

                                if didStarve > 1: larvaPsur[cohort,ind,larvaIndex,prey]=0.0
                                else:larvaPsur[cohort,ind,larvaIndex,prey]=larvaPsur[cohort,ind,larvaIndex-1,prey]*(np.exp(-mortality))
                                #
                                #print 'Prob.survival at depth : %3.2f for time %s => %3.6f for size %s'%(depth,
                                #                                                            grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                                #                                                            ,larvaPsur[cohort,ind,t,prey]*100.,
                                #                                                            L[cohort,ind,t,prey])

                                SGR[cohort,ind,larvaIndex,prey]=((W[cohort,ind,larvaIndex,prey]-W[cohort,ind,larvaIndex-1,prey])/W[cohort,ind,larvaIndex-1,prey])*100.0

                                if ind==0 and cohort==minNumberOfActiveCohort and prey==0: #releasedCohorts-1:
                                    larvaTime.append(julian/3600.)


                                julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION)
                                s1=int(startAndStop[cohort,ind,0])

                                Age[cohort,ind,prey] =julian/3600.-larvaTime[s1]
                              #  print "Age %s at time %s with SGR %s"%(Age[cohort,ind,prey],larvaIndex,  SGR[cohort,ind,larvaIndex,prey])
                                """Update time"""
                                currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                                delta=currentDate - grdSTATION.refDate

                                julian=delta.days*86400 + delta.seconds
                                t +=1
                                """This is the index that counts the time for each individual larva"""
                                larvaIndex=t-startAndStop[cohort,ind,0]

                            now=grdSTATION.refDate + datetime.timedelta(seconds=julian)
                            last=grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive)

                            startAndStop[cohort,:,1]=t

                            if now == grdSTATION.endDate or now == last:

                                if ind==(Nlarva-1) and cohort==(releasedCohorts-1) and prey==(Nprey-1):
                                    message='---> Finished IBMtime on date: %s. Writing results to file....'%(grdSTATION.refDate + datetime.timedelta(seconds=julian))
                                    showProgress(Ntime,Ntime,message)
                                    larvaTime.append(julian/3600.)

                                    IOwrite.writeStationFile(deltaH,deltaZ,grdSTATION,larvaTime,W,L,SGR,larvaTdata,larvaNauplii,larvaDepth,W_AF,larvaPsur,outputFile,startAndStop,Tlarva)
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

                            julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION)
                            hour=0
                    ind=0
            currentCohort+=1
        if loopsDone is False:
            """Show progress indicator"""
            message='---> running IBMtime-step %s of %s with %s released cohorts (%s finished)'%(t,grdSTATION.refDate + datetime.timedelta(seconds=julian),releasedCohorts,minNumberOfActiveCohort)
            showProgress(t,Ntime,message)

def real_main():

    #events = ['COLD', 'WARM']
    #events = ['REGULAR RUN']
    #events = ['CLIM RUN']
    events=  ['ESM RUN']

    if events[0]=="REGULAR RUN" or events[0]=="COLD":
        """ REGULAR run--------------------"""

        """These regular input files were created using soda2roms (extract stations)"""
        stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_NorthSea.nc",
                  "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_Iceland.nc",
                  "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_Lofoten.nc",
                  "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc"]

    if events[0]=="CLIM RUN":
        """ CLIM run--------------------"""

        """These input files were generated using soda2average/extractSODAstations.py"""
        stations = ["/Users/trond/Projects/arcwarm/SODA/soda2average/clim/climStation_NorthSea.nc",
                    "/Users/trond/Projects/arcwarm/SODA/soda2average/clim/climStation_Iceland.nc",
                    "/Users/trond/Projects/arcwarm/SODA/soda2average/clim/climStation_Lofoten.nc",
                    "/Users/trond/Projects/arcwarm/SODA/soda2average/clim/climStation_GeorgesBank.nc"]

    if events[0]=="ESM RUN":
        """ ESM run--------------------"""

        """These input files were generated using ESM2/organizeESMdata.py"""
        stations = ["/Users/trond/Projects/ESM2/ESM-northsea.nc",
                    "/Users/trond/Projects/ESM2/ESM-iceland.nc",
                    "/Users/trond/Projects/ESM2/ESM-lofoten.nc",
                    "/Users/trond/Projects/ESM2/ESM-georges.nc"]

        stations = ["/Users/trond/Projects/ESM2/ESM-northsea.nc"]


    """Notice that if you use seawifs data, the order of the stations here have to match
    The order of station seawifs data in seawifs file (see extractSeaWifs.py file)"""
    stationNames=['North Sea', 'Iceland','Lofoten', 'Georges Bank']
 

    stationNumber=0

    help.showInfo()

    print "Running ibm in mode: %s"%(events[0])
    if events[0] !='COLD':
        for station, stationName in zip(stations,stationNames):
            """REGULAR, ESM2, and CLIM RUNS ARE CALLED HERE"""
            ibm(station, stationName, stationNumber, events[0])
            stationNumber+=1

    else:
        """COLD AND WARM EVENTS"""
        for station, stationName in zip(stations,stationNames):
            for event in events:
                ibm(station, stationName, stationNumber, event)

            stationNumber+=1

def profile_main():
    # This is the main function for profiling
    # We've renamed our original main() above to real_main()
    import cProfile, pstats, StringIO, logging
    prof = cProfile.Profile()
    prof = prof.runctx("real_main()", globals(), locals())
    stream = StringIO.StringIO()
    stats = pstats.Stats(prof, stream=stream)
    stats.sort_stats("time")  # Or cumulative
    stats.print_stats(80)  # 80 = how many to print
    # The rest is optional.
    #stats.print_callees()
    #stats.print_callers()
    LOG_FILENAME = 'ibm_log.log'
    if os.path.exists(LOG_FILENAME): os.remove(LOG_FILENAME)
    logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
    logging.info("Profile data:\n%s", stream.getvalue())
    
profile_main()