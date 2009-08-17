import os, sys, string
import netCDF4
from netCDF4 import num2date, Dataset
import types
import numpy as np
import IOwrite
import datetime as datetime
"""
Import modules created as part of project
"""

""" Get self made modules"""
dir='/Users/trond/Projects/arcwarm/SODA/soda2roms'

if os.path.isdir(dir):
    sys.path.append(dir)


import grd
import IOverticalGrid
import IOtime
import IOlight
import date
import IOnetcdf
import predation

from initLarva import *

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 9)
__modified__ = datetime.datetime(2009, 8, 10)
__version__  = "1.1.5"
__status__   = "Production"

def help():
    """
    Print the helpful module docstring.
    help() -> None
    """
    
    print __doc__
    return
    
def init(station):
    """
    init: initializes the reading of files and definition of global variables
    """
    log         = True
    clim        = False
    Finish      = False
    fileNameIn  = station
    
    varlist=['temp','salt','u','v'] #,'nanophytoplankton','diatom','mesozooplankton','microzooplankton','Pzooplankton']
    startDate = datetime.datetime(1991,3,12,0,0,0)
    endDate   = datetime.datetime(1991,4,15,0,0,0)
    
    """Open the netcdf file if it existst."""
    cdf = Dataset(fileNameIn)
    
    """Calculate the sigma to meters matrix. This is important as all variables in the netcdf file are stored
    at sigma layers. To be able to convert from sigma to meters we use this function."""
    grdSTATION = grd.grdClass(fileNameIn,"STATION")

    """Get the time information and find the indices for start and stop data to extract relative to
    the time period wanted. Takes input data:
    startDate="DD/MM/YYYY" and endDate="DD/MM/YYYY" or if none given, finds all date"""
    getTimeIndices(cdf,grdSTATION,startDate,endDate,clim)
    
    """Extract the variables at the given station and store in the grdSTATION object"""
    IOnetcdf.getStationData(cdf,varlist,grdSTATION,log,clim)
    
    """Open output file:"""
    outputFile=str(station)+'_res.nc'
    if os.path.exists(outputFile): os.remove(outputFile)
    
    """Number of release dates and cohorts"""
    NReleaseDates = 1
    daysBetweenReleases=7
    """Calculate release dates for individual cohorts based on days since start date of simulations."""
    ListOfReleaseDates=[]
    for i in range(NReleaseDates):
        date=startDate+datetime.timedelta(daysBetweenReleases*(i+1))
        ListOfReleaseDates.append(date)
    
    return grdSTATION, clim, outputFile, ListOfReleaseDates, NReleaseDates, Finish

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
            depthIndex1=d; depthIndex2=d; dz1=0.5; dz2=0.5
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
    NOptDepths=abs(oldDepth-h_start) + abs(h_stop-oldDepth)
    fitnessAtDepth=0.0
    
    return oldDepth,optDepth,h_start,h_stop,NOptDepths,fitnessAtDepth

def getBehavior(stomachFullness,F,m,length,oldFitnessAtDepth,depth,optDepth):
    """Rule 4 of Behavioral Ecology paper - Kristiansen et al. 2009"""
    T=min(1.0,0.3+1000.0*(1+(length)*np.exp(length))**(-1))
    print F, m, stomachFullness, T
    if stomachFullness > T:
        beta   = 7.0
        vector = ((stomachFullness-T)/(1.0-T))**beta
    else: vector=0.0
    fitness=(((1-vector)*F) - vector*m)
    
    if fitness > oldFitnessAtDepth:
        print "Improved fitness at new depth %s vs old depth %s : new fit %s and old %s "%(optDepth,depth,fitness,oldFitnessAtDepth) 
        optDepth=depth
        oldFitnessAtDepth=fitness
                       
    return optDepth,oldFitnessAtDepth
    
def getData(julian,julianIndex,julianFileA,julianFileB,dz1,dz2,depthIndex1,depthIndex2,grdSTATION):
    """Calculate weights to use on input data from file"""
    dwA = abs(julian) - abs(julianFileA)
    dwB = abs(julianFileB) - abs(julian)
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
    
    return Tdata,Sdata,Udata,Vdata

def calculateGrowth(depth,hour,grdSTATION,dt,Wpre,Wpost,Lpost,
                    julian,julianIndex,julianFileA,julianFileB,R,enc,hand,ing,pca,Spre):
    
    """Wpost = W[cohort,ind,t,prey]
       Wpre  = W[cohort,ind,t-1,prey]
       Lpost = L[cohort,ind,t,prey]
       Lpre  = L[cohort,ind,t-1,prey]"""
       
    s_light = IOlight.surfaceLight(julian,grdSTATION.lat,hour)
    Eb = s_light*np.exp(attCoeff*(-depth))

    depthIndex1, depthIndex2, dz1, dz2 = getDepthIndex(grdSTATION,depth)
    Tdata,Sdata,Udata,Vdata = getData(julian,julianIndex,julianFileA,
                                      julianFileB,dz1,dz2,depthIndex1,
                                      depthIndex2,grdSTATION)
                        
    """Mouthsize of larvae. The larvae can only capture calanus of sizes less
    than the larval mouthsize. Folkvord et al."""
    m = np.exp (-3.27+1.818*np.log(Lpost)-0.1219*(np.log(Lpost))**2.)
    
    meta = dt*2.38e-7*np.exp(0.088*7)*((Wpre*1000)**(0.9)*0.001)
    if Eb > 0.001:
        if Lpost > 5.5: meta = (2.5*meta)
        else: meta = (1.4*meta)

    assi = 0.8*(1.0 - 0.4*np.exp(-0.002*(Wpre/mm2m-30.0)))
    GR = sec2day*np.log(0.01*(1.08 + 1.79*Tdata - 0.074*Tdata*np.log(Wpre)
                              - 0.0965*Tdata*np.log(Wpre)**2
                              + 0.0112*Tdata*np.log(Wpre)**3) + 1)*dt  
    GR_gram = (np.exp(GR) - 1)*Wpre
                   
    Em = Lpost**2/(contrast*0.1*0.2*0.75) #Size-specific sensitivity of the visual system (Fiksen,02)
    for j in range(13):
        Ap_calanus = 0.75*calanus_L1[j]*mm2m*calanus_L2[j]*mm2m
        R[j] = (IOlight.getPerceptionDistance(Em,attCoeff,Ke_larvae,Ap_calanus,Eb))*m2mm
        
    for j in range(13):
        hand[j] = 0.264*10**(7.0151*(calanus_L1[j]/Lpost)) # Walton 1992
        enc[j] = (0.667*np.pi*(R[j]**3.)*f + np.pi*(R[j]**2.)*np.sqrt(calanus_L1[j]**2.+ 2.*omega**2.)*f*tau) * calanus_D[j]* ltr2mm3
        pca[j] = enc[j]*max(0.0,min(1.0,-16.7*(calanus_L1[j]/Lpost) + 3.0/2.0))
        ing[j] = dt*pca[j]*calanus_W[j]*micro2m / (1 + hand[j]);
      
    stomachFullness =  (min(gut_size*Wpre,Spre + sum(ing)))/(Wpre*gut_size)
    
    return ing,GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness
   
def getMaxMinDepths(oldDepth,maxHourlyMove):
    if (oldDepth>maxHourlyMove and oldDepth < 100):
       h_start = oldDepth - maxHourlyMove
       h_end   = oldDepth + maxHourlyMove

    elif (oldDepth <= maxHourlyMove):
       h_start=0
       h_end=min(100,oldDepth + maxHourlyMove)
    elif (oldDepth + maxHourlyMove >= 100):
       h_start=max(1,oldDepth-maxHourlyMove)
       h_end=100
    else:
       print 'STOP: Depth error (ibm.py:getMaxMinDepths)'

    return int(h_start),int(h_end)

    
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
                    print "%s at index %i => %s"%(refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400),i,JDstart)
                    grdSTATION.startIndex=i-1
                    grdSTATION.start_date = refDate + datetime.timedelta(seconds=grdSTATION.time[i]*86400)
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
        #print 'reset time julian',julianFileA/3600.,julian/3600.,julianFileB/3600.
        julianIndex=julianIndex+1
        julianFileA=grdSTATION.time[julianIndex]*86400
        julianFileB=grdSTATION.time[julianIndex+1]*86400
       # print 'new  time julian',julianFileA/3600.,julian/3600.,julianFileB/3600.
       # print "diff to go: ",julianFileB/3600.-julian/3600.
       # print "\n"
        Finish = False
    elif julian >= grdSTATION.JDend:
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


def ibm(station):
    
    """
    Read the netcdf file and create the input arrays needed by the ibm
    """
    grdSTATION, clim, outputFile, ListOfReleaseDates, Ncohorts, Finish = init(station)
    
    """Start time in days, hours, and seconds since 1948,1,1"""
    time=grdSTATION.start_date - grdSTATION.refDate 
   
    """ julian is the time used to control the larvae, while
    julianFileA and julianFileB are the time stamps of inbetween
    the larvae recides relative to time from file. Temperature is
    interpolated to date julian from julianFileA and julianFileB"""
    #TODO: Need to also add the hours but hours is not in timedelta? No worries for soda data though
    julian=time.days*86400 + time.seconds
    julianFileA=grdSTATION.time[0]*86400
    julianFileB=grdSTATION.time[1]*86400
    Ntime=int(abs(grdSTATION.time[-1]-grdSTATION.time[0])*24.) # Number of hours to use in arrays

    if clim is True:
        julianFileA=grdSTATION.time[0]
        julianFileB=grdSTATION.time[1]
    julianIndex=0
    
    """Create array for each individual larvae"""
   
    W   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    larvaPsur   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    larvaDepth   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    W_AF=np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    S   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    Age =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    L   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
   
    ingrate  =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    Psurvive =np.ones((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    SGR      =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
    larvaTdata   =np.zeros((Ncohorts,Nlarva,Ntime,Nprey),dtype=np.float64)
   
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
    grdSTATION.Nprey=1.0
    grdSTATION.larvaPsur=1.0
    
    prey=0
    t=1
    
    print "RANDOM WEIGHT INITIALIZED"
    W[:,:,:,:]         = 0.093 
    W[:,:,:,:]         = W+ (W/4.)*np.random.random_sample(W.shape)# milligram (5mm larvae)
    W_AF[:,:,:,:]      = 0.093
    S[:,:,:,:]         = 0.3*0.06*W # 30% av max mageinnhold
    larvaTime=[]
    
    depth=5
    """==================LOOP STARTS====================="""
    """Loop over all the days of interest"""
    for day in range(Ntime):
        
        """Save the julian date and the time index so that you can reset time between each chohort and indiviual you run"""
        oldJulian=julian
        oldT=t
        """loop over all the cohorts of interest"""
        for cohort in range(Ncohorts):
            """For each day, and cohort, loop over all the individuals"""
            for ind in range(Nlarva):
                
                # TODO: Fix so that you can use any timestep. Do this by making the loop
                # count on Nsteps instead of Nhours, where Nsteps*dt=Nhours.
                
                for hour in range(Nhours):
                    
                    L[cohort,ind,t,prey] = np.exp(2.296 + 0.277*np.log(W[cohort,ind,t-1,prey]) - 0.005128*np.log(W[cohort,ind,t-1,prey])**2)
                    swimSpeed=0.261*(L[cohort,ind,t,prey]**(1.552*L[cohort,ind,t,prey]**(0.920-1.0)))-(5.289/L[cohort,ind,t,prey])
                    maxHourlyMove=(swimSpeed*(dt/1000.))/2.0 # Divided by two compared to original IBM in fortran

                    oldDepth,optDepth,h_start,h_stop,NOptDepths,oldFitnessAtDepth = initBehavior(depth,maxHourlyMove)
                   
                    for findDepth in range(NOptDepths+1):
                        depth=h_start+findDepth
                      
                        """Calculate growth at the current depth layer"""
                        ing, GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness = calculateGrowth(depth,hour,grdSTATION,dt,W[cohort,ind,t-1,prey],
                                                                                     W[cohort,ind,t,prey],L[cohort,ind,t,prey],
                                                                                     julian,julianIndex,julianFileA,julianFileB,
                                                                                     R,enc,hand,ing,pca,S[cohort,ind,t-1,prey])
                        """Calculate mortality at the current depth layer"""
                        mortality = predation.FishPred(FishDens,L[cohort,ind,t,prey]*mm2m,attCoeff,Ke_predator,Eb,dt)
                        """Calculate fitness at the current depth layer"""
                        F = sum(ing)/W[cohort,ind,t-1,prey]
                        
                        optDepth,oldFitnessAtDepth = getBehavior(stomachFullness,
                                                     F,mortality,L[cohort,ind,t,prey],
                                                     oldFitnessAtDepth,depth,optDepth)
                    
                    print "Best depth is ", optDepth, hour 
                    depth=optDepth
                    
                    """Now recalculate the growth and mortality at the optimal depth layer"""
                    ing, GR_gram,GR,meta,assi,Tdata,Eb,stomachFullness = calculateGrowth(depth,hour,grdSTATION,dt,W[cohort,ind,t-1,prey],
                                                                         W[cohort,ind,t,prey],L[cohort,ind,t,prey],
                                                                         julian,julianIndex,julianFileA,julianFileB,
                                                                         R,enc,hand,ing,pca,S[cohort,ind,t-1,prey])
                    
                    mortality = predation.FishPred(FishDens,L[cohort,ind,t,prey]*mm2m,attCoeff,Ke_predator,Eb,dt)
                    
                    S[cohort,ind,t,prey] = min(gut_size*W[cohort,ind,t-1,prey],S[cohort,ind,t-1,prey] + sum(ing[:]))

                    ingrate[cohort,ind,t,prey] = (S[cohort,ind,t,prey] - S[cohort,ind,t-1,prey])/W[cohort,ind,t-1,prey]
                    W[cohort,ind,t,prey] = W[cohort,ind,t-1,prey] + min(GR_gram + meta,S[cohort,ind,t,prey]*assi) - meta #- 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta;
                    S[cohort,ind,t,prey] = S[cohort,ind,t,prey] - ((W[cohort,ind,t,prey] - W[cohort,ind,t-1,prey]) + meta)/assi # + 0.1*act*abs(Init(l,6) - Init(l,6))/(L(l)*dt)*meta)/assi;
                    
                    W_AF[cohort,ind,t,prey] = W_AF[cohort,ind,t-1,prey] + (np.exp(GR)-1)*W_AF[cohort,ind,t-1,prey]
                    larvaTdata[cohort,ind,t,prey] = Tdata
                    larvaDepth[cohort,ind,t,prey] = depth
                    #Psurvive[ind,hour] = np.exp(-a*(L[ind]**b) - C2*(1-Pe)*Rpisci**2)
                    #for s in range(Nhours):
                    #TODO: FIX the survival probability
                    larvaPsur[cohort,ind,t,prey]=larvaPsur[cohort,ind,t-1,prey]*(np.exp(-Psurvive[cohort,ind,t,prey]))
                    #print 'Probability of survival at depth : %3.2f => %6.2f'%(depth,grdSTATION.larvaPsur*100.)
                    
                    SGR[cohort,ind,t,prey]=((W[cohort,ind,t,prey]-W[cohort,ind,t-1,prey])/W[cohort,ind,t-1,prey])*100.0
                        
                    """Update time"""
                    if ind==(Nlarva-1) and cohort==(Ncohorts-1): 
                        julianFileA, julianFileB,julianIndex, Finish = updateFileIndices(julian,julianIndex,julianFileB,julianFileA,grdSTATION,clim)
                        larvaTime.append(julian/3600.)
            
                       # print "Added time to list", julian/3600.
                        
                    currentDate=grdSTATION.refDate + datetime.timedelta(seconds=julian + dt)
                    delta=currentDate - grdSTATION.refDate 
                    
                    julian=delta.days*86400 + delta.seconds
                    t +=1
                
                 
                # TODO: FIX these lines if you want to use climatological station data 
                #if depth <= Nlarva-1:
                #    
                #    if clim is False:   
                #        print_date=str(currentDate[0])+'-'+print_mm+'-'+str(print_dd)+'T%2i:00'%(hour)
                #    else:
                #        print_date='CLIM-'+print_mm+'-'+str(print_dd)+'T%2i:00'%(hour)
                
                """Calculate 24 hour SGR """        
                max_g = 1.08 + 1.79 * Tdata - 0.074 * Tdata*np.log(W[cohort,ind,t,prey]) - 0.0965 * Tdata *((np.log(W[cohort,ind,t,prey]))**2) + 0.0112 * Tdata *((np.log(W[cohort,ind,t,prey])**3.))
                     
                if ind==(Nlarva-1) and cohort==(Ncohorts-1): 
                    if Finish is True:
                
                        """Write the results to file"""
                        IOwrite.writeStationFile(grdSTATION,larvaTime,W,L,SGR,larvaTdata,larvaDepth,W_AF,aveLight,larvaPsur,outputFile)
                        os.sys.exit()
                        
                """ Resetting values for all larvae"""
                #if ind==(Nlarva-1): 
                #    for i in range(Nlarva):
                #        W[cohort,ind,t,prey]    = 0.093 # milligram (5mm larvae)
                #        W_AF[cohort,ind,t,prey] = 0.093
                #        S[cohort,ind,t,prey]    = 0.3*0.06*W[cohort,ind,t,prey] # 30% av max mageinnhold
                #        Age[cohort,ind,t,prey] = 0
                #        Psurvive[cohort,ind,t,prey]=1.0
                #        grdSTATION.larvaPsur=1.0
                if ind!=(Nlarva-1) and cohort!=Ncohorts:       
                    julian=oldJulian
                    t=oldT
                hour=0
        """Show progress indicator"""
        message='---> running IBM'
        p=int( ((t*1.0+2)/(1.0*Ntime))*100.)
        progress.render(p,message)
                        
            #print t, cohort, ind, currentDate, W[cohort,ind,t,prey]
    
    
if __name__=="__main__":
   
    try:
        import psyco
        psyco.log()
        psyco.full(memory=100)
        psyco.profile(0.05, memory=100)
        psyco.profile(0.2)
    except ImportError:
        pass
    
    stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_EastandWestGreenland.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_lofoten.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_NorthSea.nc",
              "/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_Iceland.nc"]
    stations=["/Users/trond/Projects/arcwarm/SODA/soda2roms/stations/station_GeorgesBank.nc"]

    
    for station in stations:
        
        ibm(station)
  