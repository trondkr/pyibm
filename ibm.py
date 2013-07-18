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
from initLarva import *
"""Import f2py Fortran modules"""
import calclight
import perception
import bioenergetics
import predF90
import IOstation
from memoize import memoize
import IOdata
import behavior

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 9)
__modified__ = datetime.datetime(2010, 2, 10)
__version__  = "1.1.5"
__status__   = "Production"
__revisions__= "8.6.09, 2.12.09, 19.01.10, 26.01.10, 10.02.2010"

def revisions():
    """
    03.10.2011 TK Changed the number of prey species to 12 to speedup the calculations
                  Found a bug that caused zooplankton to not be calculated sometimes when
                  event was set to "None". Unclear how important that bug was.
    """

def infoOnResolution(grdSTATION):

    print 'Defined resolution of simulation:'
    print 'Temporal : %s hours per time-steps (%s hours)'%(deltaH, grdSTATION.refDate +
                                                             datetime.timedelta(seconds=deltaH*3600.)
                                                             - grdSTATION.refDate)
    print 'Vertical : %s cm (rounds to %i decimalpoints)'%(deltaZ*100.,lastDecimal)
    print '-----------------------------------------------------------------\n'

def infoOnOptions(grdSTATION):
    if "seawifs" in grdSTATION.OPTIONS:
        doc="""\nSEAWIFS IS ON: The model is run using seawifs monthly values
        as input to the model to calculate zooplankton values. This is the
        standard for the COLD/WARM simulations where chlorophyll-a values
        are converted into zooplankton.\n"""
        print doc

    if "foodLimited" in grdSTATION.OPTIONS:
        doc="""\nUSEFOODLIMITED IS ON: This is the default option for all
        runs as it assumes that the larval growth rates are food limited. If you
        turn off this option (remove from OPTIONS list), the growth rates are
        calculated assuming unlimited growth rates according to Folkvord 2005.
        This is the deafult for the Tristan Rouyer simulations.\n"""
    else:
        doc="""\nSEFOODLIMITED IS OFF: This means that growth rates are
        estimated assuming unlimited food supply for the larvae. This calculates
        growth rates according to Folkvord 2005.\n"""
        print doc

    if "useAverageFile" in grdSTATION.OPTIONS:
        doc="""\nSUSEAVERAGEFILE IS ON: The model is run using seawifs
        together with temperature anomalies monthly values to estimate
        zooplankton values as input to the model. This is the
        standard for the COLD/WARM simulations where chlorophyll-a values and
        temperature anomalies are converted into zooplankton. This option requires
        the climatology file:
        "/Users/trond/Projects/arcwarm/SODA/soda2average/clim/averageSODA1961-1990.nc"
        to exists.\n"""
        print doc
    if "useZooplanktonModel" in grdSTATION.OPTIONS:
        doc="""Calculate the zooplankton densities using large and small phytoplankton as input.
        This approach assumes that the order of the variables is:
        varlist=['temp','salt','nh4sm','nh4lg','no3sm','no3lg','chla','taux','tauy']
        If this is not correct, change the order in init function and re-run.
        """
        print doc
    if "useDetrendedTemperature" in grdSTATION.OPTIONS:
        doc="""The predicted temperature trend has been removed from the temperature
        data used as input to the IBM. This is a detrended run used to estimate the
        relative effect of having both predicted changes in temperature and phytoplankton
        versus only predicted changes in phytoplankton and using detrended temperature. The
        temperature is detrended using the Python script that creates the input files to the IBM
        form the ESM outputs: /Users/trond/Projects/ESM2/organize5DEGESMdata.py\n"""
        print doc

def init(station,stationName,event,eventOption,historyRun):
    """
    init: initializes the reading of files and definition of global variables
    """
    log            = True
    Finish         = False
    fileNameIn     = station

    if event=="ESM RUN":
      #  OPTIONS=["foodLimited", "useZooplanktonModel","biasCorrectTemperature","useDetrendedTemperature","useReducedFishDensity"]
        OPTIONS=["foodLimited", "useZooplanktonModel","biasCorrectTemperature"]

    """Option for using seawifs data as prey abundance"""
    seawifsFileName="/Users/trond/Projects/seawifs/chlo-stations.nc"
    if "seawifs" in OPTIONS:
        seawifsArray=IOnetcdf.getSeaWifs(seawifsFileName)
    else:
        seawifsArray=None

    """Open the netcdf file if it existst."""
    if os.path.exists(fileNameIn):
        cdf = Dataset(fileNameIn)
    else: print "File %s does not exist"%(fileNameIn)

    """Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function."""
    grdSTATION = grd.grdClass(fileNameIn,"STATION")

    if event=="ESM RUN":
        if "useDetrendedTemperature" in OPTIONS:
            print "NOTE: RUNNING WITH DETRENDED TEMPERATURE!"
            varlist=['detemp','salt','nh4sm','nh4lg','no3sm','no3lg','chla','taux','tauy']
        else:
            varlist=['temp','salt','nh4sm','nh4lg','no3sm','no3lg','chla','taux','tauy']
        if historyRun is False:
            startDate = datetime.datetime(2001,1,15,1,0,0)
            endDate   = datetime.datetime(2099,12,31,1,0,0)
        else:
            startDate = datetime.datetime(1950,1,1,1,0,0)
            endDate   = datetime.datetime(1998,12,31,1,0,0)
        grdSTATION.chlaValue=0.0

    """Define the list of options used during excecution"""
    grdSTATION.OPTIONS=OPTIONS
    grdSTATION.log=log

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
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,stationName,eventOption)
    """Store the maximum and minimum temperature values for the time period of interest, and
    use these to find the prey density as a function of temperature: function calculateGrowth"""

    grdSTATION.minT = grdSTATION.data[:,:,0].min()
    grdSTATION.maxT = grdSTATION.data[:,:,0].max()

    """Open output file:"""
    outputFile="IBM_"+os.path.basename(station)

    if "useDetrendedTemperature" in OPTIONS:
        outputFile="IBM_deTemp_"+os.path.basename(station)

    if "useReducedFishDensity" in grdSTATION.OPTIONS:
        if not os.path.exists("resultsReducedFishDens"): os.mkdir("resultsReducedFishDens")
        if os.path.exists(outputFile): os.remove(outputFile)
        if os.path.exists("resultsReducedFishDens/"+outputFile): os.remove("resultsReducedFishDens/"+outputFile)
    else:
        if not os.path.exists("results"): os.mkdir("results")
        if os.path.exists(outputFile): os.remove(outputFile)
        if os.path.exists("results/"+outputFile): os.remove("results/"+outputFile)


    print "Outputfile will be named: %s"%(outputFile)

    if not os.path.exists("resultsReducedFishDens"): os.mkdir("resultsReducedFishDens")
    if os.path.exists(outputFile): os.remove(outputFile)
    if os.path.exists("resultsReducedFishDens/"+outputFile): os.remove("resultsReducedFishDens/"+outputFile)

    listOfReleaseDates=[]
    """Calculate release dates for individual cohorts based on days since start date of simulations."""
    listOfReleaseDates.append(startDate)

    for i in xrange(NReleaseDatesInit):
        date=startDate+datetime.timedelta(daysBetweenReleases*(i+1))
        if endDate >= date + datetime.timedelta(days=NDaysAlive):
            listOfReleaseDates.append(date)

    NReleaseDates=len(listOfReleaseDates)
    grdSTATION.listOfReleaseDates=listOfReleaseDates

    print "\nThis simulation will release a total of %s cohorts between the following dates:"%(NReleaseDates)
    #for h in xrange(NReleaseDates):
    print "--> %s"%(listOfReleaseDates[0])
    print "--> %s"%(listOfReleaseDates[-1])
    print "\n"

    print "IBM will run from date %s to %s\n"%(grdSTATION.listOfReleaseDates[0] + datetime.timedelta(days=0),
                                              grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive))
    print "-----------------------------------------------------------------"

    grdSTATION.stationNumber=0

    infoOnResolution(grdSTATION)
    infoOnOptions(grdSTATION)

    return grdSTATION, outputFile, NReleaseDates, Finish, seawifsArray

def findActiveCohorts(isDead,isReleased):
    co=sum(isDead)
    co2=sum(isReleased)
    activeCohorts=abs(co2-co)
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
    depthArray=np.abs(grdSTATION.depth)
    depth=abs(depth)

    while d < len(depthArray) and depthFOUND==False:

        if depthArray[d] <= depth <= depthArray[d+1] and depthFOUND==False:
            dz1=1.0-abs(depthArray[d]-depth)/abs(depthArray[d]-depthArray[d+1])
            dz2=1.0-abs(depthArray[d+1]-depth)/abs(depthArray[d]-depthArray[d+1])
            depthIndex1=d; depthIndex2=d+1
            depthFOUND=True
            #print 'depth 1:',depthArray[d], depth, depthArray[d+1],depthIndex1, depthIndex2

        elif depthArray[0] > depth and depthFOUND==False:
            depthIndex1=0; depthIndex2=0; dz1=0.5; dz2=0.5
            depthFOUND=True
            #print 'depth 2:',grdSTATION.depth[d], depth,depthIndex1, depthIndex2

        elif depthArray[-1] < depth and depthFOUND==False:
            depthIndex1=-1; depthIndex2=-1; dz1=0.5; dz2=0.5
            depthFOUND=True
            #print 'depth 3:',grdSTATION.depth[-1], depth,depthIndex1, depthIndex2
        elif d==len(depthArray)-1 and depthFOUND==False:
            print "No valid index for depth (%sm) found"%(depth)

        d+=1

    return depthIndex1, depthIndex2, dz1, dz2

def initBehavior(depth,maxHourlyMove,grdSTATION):
    oldDepth=depth  #The previous depth for a given larvae for a time step inside optimal loop
    optDepth=depth  # The optimal depth for a given larvae for a time step
    h_start, h_stop = getMaxMinDepths(oldDepth,maxHourlyMove,grdSTATION.deepestDepthAllowed)

    NOptDepths=int((abs(oldDepth-h_start) + abs(h_stop-oldDepth) )/float(deltaZ))
    # Limit the number of depths to check to 6 (remove this will reduce speed)
    NOptDepths=min(10,NOptDepths)
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

def calculateLength(Larval_wgt, Larval_m):
    oldLarval_m = Larval_m
    Larval_m = np.exp(2.296 + 0.277*np.log(Larval_wgt) - 0.005128*np.log(Larval_wgt)**2)
    if Larval_m < oldLarval_m:
        #print 'Kept the old length', oldLarval_m, Larval_m
        Larval_m = oldLarval_m

    return Larval_m

def getMaxMinDepths(oldDepth,maxHourlyMove,deepestDepthAllowed):
    """The value you give deepstDepthAllowed will be the absolute deepest
    depth the larvae will move to. Useful to set this depth equal to bottom.
    grdSTATION.deepestDepthAllowed is defined from data in function IOnetcdf.getData"""

    if (oldDepth>maxHourlyMove and oldDepth+maxHourlyMove < deepestDepthAllowed):
        h_start = oldDepth - maxHourlyMove
        h_end   = oldDepth + maxHourlyMove

    elif (oldDepth <= maxHourlyMove):
        h_start=0.0
        h_end=min(deepestDepthAllowed,oldDepth + maxHourlyMove)

    elif (oldDepth + maxHourlyMove >= deepestDepthAllowed):
        h_start=max(1,oldDepth-maxHourlyMove)
        h_end=deepestDepthAllowed
    else:
        print "getMaxMinDepths",oldDepth,maxHourlyMove
        print 'STOP: Depth error (ibm.py:getMaxMinDepths)'

    return h_start, h_end

def calculateAttenuationCoeff(Chla):
    """See Fiksen et al. 2002 for details (relationship taken from Riley 1956)
    units of Chla is mg/m3"""
    k0 = 0.1
    return k0 + 0.054*Chla**(2./3.) + 0.0088*Chla

def getTimeIndices(cdf,grdSTATION,startDate,endDate,event):
    """Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """

    if event=="ESM RUN":
        if startDate.year < 2000:
            refDate    = datetime.datetime(1861,1,1,0,0,0)
        else:
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


def updateFileIndices(julian,julianIndex,julianFileB,julianFileA,myTime,myJDend,startDate,event):
    """This routine makes sure that we are reading at the correct time indices
    in the input data. If we have moved forward in time enough to move past the current time window of two
    time stamps, then update the indices and give the new indices.
    Remember that input time format from netCDF files must be
    'Days since 1948.1.1'
    """

    if event=="ESM RUN":
        if startDate.year < 2000:
            refDate    = datetime.datetime(1861,1,1,0,0,0)
        else:
            refDate    = datetime.datetime(2001,1,1,0,0,0)

    #print "searching: \n%s \n%s \n%s"%(refDate + datetime.timedelta(seconds=julianFileA),
    #                                   refDate + datetime.timedelta(seconds=julian),
    #                                   refDate + datetime.timedelta(seconds=julianFileB))

    if julian >=julianFileB and julian < myJDend:
        julianIndex=julianIndex+1
        julianFileA=myTime[julianIndex]*86400
        julianFileB=myTime[julianIndex+1]*86400
        Finish = False
    elif julian < julianFileB and julian < julianFileA:
        """Now, if you run with more than one individual you have to reset the
        fileIndices for where in the input files you interpolate from between each individual.
        Therefore, here we search backwards to find the correct start and stop indices for the new
        individual."""

        while julian < julianFileB and julian < julianFileA:
            julianIndex=julianIndex-1
            julianFileA=myTime[julianIndex]*86400
            julianFileB=myTime[julianIndex+1]*86400
            Finish = False

    elif julian > myJDend:
        Finish = True
    else:
        Finish = False
    julianFileA=myTime[julianIndex]*86400
    julianFileB=myTime[julianIndex+1]*86400

    return julianFileA, julianFileB, julianIndex, Finish


def showProgress(t,Ntime,message):
    p=int( ((t*1.0)/(1.0*(Ntime-24*dt_per_day)))*100.)
    progress.render(p,message)

def ibm(station,stationName,stationNumber,event,eventOption,historyRun):

    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION, outputFile, Ncohorts, Finish, seawifsArray = init(station,stationName,event,eventOption,historyRun)

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
    W            =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    larvaPsur    =(np.ones((Ncohorts,Nlarva,Tlarva),dtype=np.float64))
    larvaDepth   =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    W_AF         =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    S            =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    Age          =(np.ones((Ncohorts,Nlarva),dtype=np.float64))*missingValue
    startAndStop =(np.ones((Ncohorts,Nlarva,2),dtype=np.float64))*missingValue
    L            =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue

    ingrate      =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue

    """Store mortality trates from fish, invertebrates, and occurrence of starvation"""
    Mfish        =(np.zeros((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    Minve        =(np.zeros((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    Mstar        =(np.zeros((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue

    SGR          =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    larvaTdata   =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue
    larvaNauplii =(np.empty((Ncohorts,Nlarva,Tlarva),dtype=np.float64))*missingValue

    R=np.empty(12,dtype=np.float64)
    enc=np.empty(12,dtype=np.float64)
    hand=np.empty(12,dtype=np.float64)
    ing=np.empty(12,dtype=np.float64)

    grdSTATION.Initialized=False
    grdSTATION.Ndays=Ntime/24.
    grdSTATION.Nlarva=Nlarva
    grdSTATION.Ncohorts=Ncohorts
    grdSTATION.larvaPsur=1.0
    grdSTATION.firstBatch=True
    grdSTATION.firstRun=True

    t=1
    loopsDone = False
    grdSTATION.dayOfYear = -9
    if randomWgt==1: print "RANDOM WEIGHT INITIALIZED"

    W[:,:,:]         = initWgt + ((initWgt*0.3)* (np.random.random_sample(W.shape))*np.random.random_integers(-1,1,W.shape)*randomWgt) # milligram (5mm larvae)
    W_AF[:,:,:]      = initWgt
    S[:,:,:]         = stomach_threshold*gut_size*W # 30% av max mageinnhold
    L[:,:,:]         = calculateLength(initWgt, 0.0)

    larvaTime=[]

    depth=initDepth
    localFishDens=FishDens
    savedYear=-9
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

        if "seawifs" in grdSTATION.OPTIONS:
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

                for ind in xrange(Nlarva):
                    if day==0:
                        depth =initDepth * np.random.random_sample()
                        print "initdepth is %s and depth is %s"%(initDepth,depth)
                    if loopsDone : break
                    """We create an array that controls the entrance and exit points in time for
                    a given larvae"""

                    if startAndStop[cohort,ind,0]<=missingValue:
                        startAndStop[cohort,ind,0]=t-1
                        Age[cohort,ind]=0
                    """This is the index that counts the time for each individual larva"""
                    larvaIndex=int(t-startAndStop[cohort,ind,0])

                    isAliveBool = isAlive(julian,Age[cohort,ind],cohort,isDead)
                    if isAliveBool == False:

                        if noOneLeft==True and ind==(Nlarva-1) and cohort==(releasedCohorts-1):
                            currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*24.)
                            delta=currentDate - grdSTATION.refDate
                            julian=delta.days*86400 + delta.seconds
                            noOneLeft=True
                       #     print "Updating time as no cohort is alive currently",currentDate
                    else:
                        noOneLeft=False
                        """We calculate swimspeed once every 24 hours and assume it does not change during that time-period"""
                        swimSpeed=0.261*(L[cohort,ind,larvaIndex-1]**(1.552*L[cohort,ind,larvaIndex-1]**(0.920-1.0)))-(5.289/L[cohort,ind,larvaIndex-1])
                       # if np.isnan(swimSpeed):
                       #     print "isnan",cohort,ind,larvaIndex,t, L[cohort,ind,larvaIndex-1],L[cohort,ind,larvaIndex-2],L[cohort,ind,larvaIndex-3]
                       #     print "W,L",W[cohort,ind,larvaIndex-1],W[cohort,ind,larvaIndex-2]
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

                            if event=="ESM RUN":
                                attCoeff=calculateAttenuationCoeff(grdSTATION.chlaValue)
                                beamAttCoeff=attCoeff*3.0

                            """OPTIMAL DEPTH == Find the optimal depth to stay at given ratio of ingestion and mortality rate"""
                            for depth in sampleDepths:

                                diffZ=abs(oldDepth - depth)

                                """Light level at current depth:"""
                                Eb = surfaceLight*np.exp(attCoeff*(-depth))
                                """Calculate growth and ingestion at the current depth:"""

                                if event=="ESM RUN":
                                    zoop=np.empty(len(grdSTATION.zooplankton[0,0,:]))

                                # TODO: NOTE THAT grdSTATION.zooplankton is currently not
                                # defined for any other run than ESM RUN!

                                """==>INDEX calculations==============================="""
                                depthIndex1, depthIndex2, dz1, dz2 = getDepthIndex(grdSTATION,depth)

                                """==> EXTRACT data from files and interpolate to depth and time"""
                                Tdata=0.0;Sdata=0.0;NH4SMdata=0.0;NH4LGdata=0.0;NO3SMdata=0.0;NO3LGdata=0.0;CHLAdata=0.0;surfaceCHLA=0.0;windX=0.0;windY=0.0
                                #print ind, hour, depth, julian,julianIndex+1,julianFileA+1,julianFileB+1
                                Tdata,Sdata,NH4SMdata,NH4LGdata,NO3SMdata,NO3LGdata,CHLAdata,surfaceCHLA,windX,windY,zoop=IOdata.io.getdata(Tdata,Sdata,NH4SMdata,
                                                                                        NH4LGdata,NO3SMdata,NO3LGdata,CHLAdata,surfaceCHLA,windX,
                                                                                        windY,zoop,julian,julianIndex+1,julianFileA+1,
                                                                                        julianFileB+1,dz1,dz2,depthIndex1+1,depthIndex2+1,
                                                                                        grdSTATION.data,
                                                                                        grdSTATION.dataXY,
                                                                                        grdSTATION.zooplankton,event,
                                                                                        len(grdSTATION.data[:,0,0]),len(grdSTATION.data[0,:,0]),
                                                                                        len(grdSTATION.data[0,0,:]),len(grdSTATION.dataXY[0,:]),
                                                                                        len(zoop))


                                #NOTE: The chlorophyll is now taken from surface
                                grdSTATION.chlaValue=surfaceCHLA

                                """==> GROWTH RATE calculations at the given depth and time"""
                                suming=0.0; ing[:]=0; GR_mg=0.0; 0.0; GR=0.0; meta=0.0; assi=0.0; stomachFullness=0.0

                                suming,ing,GR_mg,meta,assi,stomachFullness = bioenergetics.bioenergetics.growth(suming,ing,GR_mg,meta,assi,stomachFullness,
                                                                                zoop,Tdata,windX,windY,S[cohort,ind,larvaIndex-1],prey_AREA,prey_LENGTH,
                                                                                prey_D,prey_WGT,L[cohort,ind,larvaIndex-1],W[cohort,ind,larvaIndex-1],
                                                                                sec2day, mg2ug,ltr2mm3, m2mm,mm2m,f,tau,ug2mg,Eb,contrast,Ke_larvae,
                                                                                beamAttCoeff,MultiplyPrey,deltaH,dt,depth,gut_size, len(zoop))

                                """Calculate mortality at the current depth layer"""
                                myMortality=0.0; Starved=0.0; StarvationMortality=0.0; FishMortality=0.0; InvertebrateMortality=0.0
                                mortality, didStarve,StarvationMortality, FishMortality, InvertebrateMortality = predF90.predation.fishpredandstarvation(myMortality,Starved,
                                                                                            StarvationMortality, FishMortality, InvertebrateMortality,
                                                                                            localFishDens,L[cohort,ind,larvaIndex-1]*mm2m,
                                                                                            W[cohort,ind,larvaIndex-1],
                                                                                            Eb,suming,stomachFullness,beamAttCoeff,m2mm,
                                                                                            Ke_predator,deadThreshold,deltaH,dt)

                                #print "Depth: %s Temp: %s mortality: %s suming: %s"%(depth,Tdata,mortality,suming)

                                """Calculate the cost of being active. The more you swim the more you use energy in raltive ratio to
                                routine metabolims. Trond Kristiansen, 23.03.2010"""
                                activityCost= (diffZ/maxDiffZ)*(meta*costRateOfMetabolism)


                                """Calculate fitness at the current depth layer"""
                                F = (suming-activityCost)/W[cohort,ind,larvaIndex-1]

                                optDepth,oldFitness = behavior.behavior.getbehavior(stomachFullness,
                                                             F,mortality,L[cohort,ind,larvaIndex-1],
                                                             depth,optDepth,oldFitness)


                            """Set the optimal depth equal to result of optimal loop and update light at depth"""
                            diffZ=abs(depth-optDepth)
                            depth=optDepth

                            Eb = surfaceLight*np.exp(attCoeff*(-depth))
                            """Now recalculate the growth and mortality at the optimal depth layer"""

                            """Extract the environmental data for the optimal depth level"""
                            depthIndex1, depthIndex2, dz1, dz2 = getDepthIndex(grdSTATION,depth)

                            """==> EXTRACT data from files and interpolate to depth and time"""
                            Tdata=0.0;Sdata=0.0;NH4SMdata=0.0;NH4LGdata=0.0;NO3SMdata=0.0;NO3LGdata=0.0;CHLAdata=0.0;surfaceCHLA=0.0;windX=0.0;windY=0.0

                            Tdata,Sdata,NH4SMdata,NH4LGdata,NO3SMdata,NO3LGdata,CHLAdata,surfaceCHLA,windX,windY,zoop=IOdata.io.getdata(Tdata,Sdata,NH4SMdata,
                                                                                    NH4LGdata,NO3SMdata,NO3LGdata,CHLAdata,surfaceCHLA,windX,
                                                                                    windY,zoop,julian,julianIndex+1,julianFileA+1,
                                                                                    julianFileB+1,dz1,dz2,depthIndex1+1,depthIndex2+1,
                                                                                    grdSTATION.data,
                                                                                    grdSTATION.dataXY,
                                                                                    grdSTATION.zooplankton,event,
                                                                                    len(grdSTATION.data[:,0,0]),len(grdSTATION.data[0,:,0]),
                                                                                    len(grdSTATION.data[0,0,:]),len(grdSTATION.dataXY[0,:]),
                                                                                    len(zoop))
                            grdSTATION.chlaValue=surfaceCHLA

                            """==> GROWTH RATE calculations at the given depth and time"""
                            suming=0.0; ing[:]=0; GR_mg=0.0; 0.0; GR=0.0; meta=0.0; assi=0.0; stomachFullness=0.0

                            suming,ing,GR_mg,meta,assi,stomachFullness = bioenergetics.bioenergetics.growth(suming,ing,GR_mg,meta,assi,stomachFullness,
                                                                                zoop,Tdata,windX,windY,S[cohort,ind,larvaIndex-1],prey_AREA,prey_LENGTH,
                                                                                prey_D,prey_WGT,L[cohort,ind,larvaIndex-1],W[cohort,ind,larvaIndex-1],
                                                                                sec2day, mg2ug,ltr2mm3, m2mm,mm2m,f,tau,ug2mg,Eb,contrast,Ke_larvae,
                                                                                beamAttCoeff,MultiplyPrey,deltaH,dt,depth,gut_size, len(zoop))



                            myMortality=0.0; Starved=0.0; StarvationMortality=0.0; FishMortality=0.0; InvertebrateMortality=0.0
                            mortality, didStarve,StarvationMortality, FishMortality, InvertebrateMortality = predF90.predation.fishpredandstarvation(myMortality,Starved,
                                                                                            StarvationMortality, FishMortality, InvertebrateMortality,
                                                                                            localFishDens,L[cohort,ind,larvaIndex-1]*mm2m,
                                                                                            W[cohort,ind,larvaIndex-1],
                                                                                            Eb,suming,stomachFullness,beamAttCoeff,m2mm,
                                                                                            Ke_predator,deadThreshold,deltaH,dt)

                            """Calculate the cost of the activity performed during last hour for vertical behavior"""
                            activityCost= (diffZ/maxDiffZ)*(meta*costRateOfMetabolism)

                            S[cohort,ind,larvaIndex] = max(0.0, min(gut_size*W[cohort,ind,larvaIndex-1],S[cohort,ind,larvaIndex-1] + sum(ing[:])))
                            ingrate[cohort,ind,larvaIndex] = max(0.0, (S[cohort,ind,larvaIndex] - S[cohort,ind,larvaIndex-1])/W[cohort,ind,larvaIndex-1])

                            if "foodLimited" in grdSTATION.OPTIONS:
                                W[cohort,ind,larvaIndex] = W[cohort,ind,larvaIndex-1] + min(GR_mg + meta,S[cohort,ind,larvaIndex]*assi) - meta - activityCost
                            else:
                                """Run as food unlimited: assuming growth is according to Folkvord 2005"""
                                W[cohort,ind,larvaIndex] = W[cohort,ind,larvaIndex-1] + GR_mg - activityCost

                            S[cohort,ind,larvaIndex] = max(0.0, S[cohort,ind,larvaIndex] - ((W[cohort,ind,larvaIndex] - W[cohort,ind,larvaIndex-1]) + meta)/assi)
                            L[cohort,ind,larvaIndex] = calculateLength(W[cohort,ind,larvaIndex], L[cohort,ind,larvaIndex-1])
                            W_AF[cohort,ind,larvaIndex] = W_AF[cohort,ind,larvaIndex-1] + (np.exp(GR)-1)*W_AF[cohort,ind,larvaIndex-1]

                            larvaTdata[cohort,ind,larvaIndex] = Tdata
                            larvaDepth[cohort,ind,larvaIndex] = depth
                            larvaNauplii[cohort,ind,larvaIndex] = sum(zoop)

                            Mfish[cohort,ind,larvaIndex]=FishMortality
                            Minve[cohort,ind,larvaIndex]=InvertebrateMortality
                            Mstar[cohort,ind,larvaIndex]=didStarve

                            #print grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH), Tdata, depth

                            if didStarve > 1000: larvaPsur[cohort,ind,larvaIndex]=0.0
                            else:larvaPsur[cohort,ind,larvaIndex]=larvaPsur[cohort,ind,larvaIndex-1]*(np.exp(-mortality))

                            #print 'Prob.survival %s at depth : %3.2f for time %s => %3.6f for size %s'%(cohort,depth,
                            #                                                            grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                            #                                                            ,larvaPsur[cohort,ind,larvaIndex]*100.,
                            #                                                            L[cohort,ind,larvaIndex])

                            SGR[cohort,ind,larvaIndex]=((W[cohort,ind,larvaIndex]-W[cohort,ind,larvaIndex-1])/W[cohort,ind,larvaIndex-1])*100.0

                            if ind==0 and cohort==minNumberOfActiveCohort:
                                larvaTime.append(julian/3600.)

                            julianFileA,julianFileB,julianIndex,Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION.time,
                                                                                           grdSTATION.JDend,grdSTATION.startDate,event)
                            s1=int(startAndStop[cohort,ind,0])

                            Age[cohort,ind] =julian/3600.-larvaTime[s1]
                          #  print "Age %s at time %s with SGR %s"%(Age[cohort,ind,prey],larvaIndex,  SGR[cohort,ind,larvaIndex,prey])
                            """Update time"""
                            currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt*deltaH)
                            delta=currentDate - grdSTATION.refDate

                            if currentDate.year > savedYear:
                                if "useReducedFishDensity" in grdSTATION.OPTIONS:
                                    localFishDens=(localFishDens-2.4e-7) #16% reduction over 100 years
                                    scaleInvertebrate=scaleInvertebrate+2.4e-7
                                    print "Reduced invertebrate: %s and fishdensity to: %s (years %s to %s)"%(scaleInvertebrate,localFishDens,savedYear,currentDate.year)
                                    savedYear=currentDate.year

                            julian=delta.days*86400 + delta.seconds
                            t +=1
                            """This is the index that counts the time for each individual larva"""
                            larvaIndex=int(t-startAndStop[cohort,ind,0])

                        now=grdSTATION.refDate + datetime.timedelta(seconds=julian)
                        last=grdSTATION.listOfReleaseDates[-1] + datetime.timedelta(days=NDaysAlive)

                        startAndStop[cohort,:,1]=t

                        if now == grdSTATION.endDate or now == last:

                            if ind==(Nlarva-1) and cohort==(releasedCohorts-1):
                                message='---> Finished IBMtime on date: %s. Writing results to file....'%(grdSTATION.refDate + datetime.timedelta(seconds=julian))
                                showProgress(Ntime,Ntime,message)
                                larvaTime.append(julian/3600.)

                                IOwrite.writeStationFile(deltaH,deltaZ,grdSTATION,larvaTime,W,L,SGR,larvaTdata,larvaNauplii,larvaDepth,W_AF,larvaPsur,outputFile,startAndStop,Tlarva,Mfish,Minve,Mstar)

                                if "useReducedFishDensity" in grdSTATION.OPTIONS:
                                    move(outputFile, "resultsReducedFishDens/"+outputFile)
                                else:
                                    move(outputFile, "results/"+outputFile)
                                loopsDone = True

                                print "Result file is saved as %s"%("resultsReducedFishDens/"+outputFile)
                                break

                        if ind==(Nlarva-1) and cohort==(releasedCohorts-1):
                            julian=julian
                            t=t
                        else:
                            julian=oldJulian
                            t=oldT

                        julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION.time,
                                                                                         grdSTATION.JDend,grdSTATION.startDate,event)
                        hour=0
                ind=0
            currentCohort+=1
        if loopsDone is False:
            """Show progress indicator"""
            message='---> running IBMtime-step %s of %s with %s released cohorts (%s finished)'%(t,grdSTATION.refDate + datetime.timedelta(seconds=julian),releasedCohorts,minNumberOfActiveCohort)
            showProgress(t,Ntime,message)

def real_main():

    events  = ['ESM RUN']

    for historyRun in [False]:
        if events[0]=="ESM RUN":
            """ ESM run--------------------"""
            if historyRun is False:
                """These input files were generated using ESM2/organize5DEGESMdata.py"""
                stations = ["/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_northsea.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_iceland.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_lofoten.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_georges.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_wglnd.nc"]

            #    stations=["/Users/trond/Projects/ESM2/organized5deg/ESM_2001_2100_georges.nc"]

            else:
                stations = ["/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_northsea.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_iceland.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_lofoten.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_georges.nc",
                            "/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_wglnd.nc"]

           #     stations = ["/Users/trond/Projects/ESM2/organized5deg/ESM_1941_2000_georges.nc"]


        """Notice that if you use seawifs data, the order of the stations here have to match
        The order of station seawifs data in seawifs file (see extractSeaWifs.py file)"""
        stationNames=['North Sea', 'Iceland','Lofoten', 'Georges Bank', 'West Greenland']
      #  stationNames=['Georges Bank']

        stationNumber=0

        help.showInfo()

        print "Running ibm in mode: %s"%(events[0])
        if events[0] in ['REGULAR RUN','CLIM RUN','ESM RUN']:
            eventOption=None
            for station, stationName in zip(stations,stationNames):
                """REGULAR, ESM2, and CLIM RUNS ARE CALLED HERE"""
                ibm(station, stationName, stationNumber, events[0], eventOption,historyRun)
                stationNumber+=1


def profile_main():
    PROFILE=False
    if PROFILE is True:
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
    else:
        real_main()

profile_main()
