#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2010, 1, 26)
__modified__ = datetime.datetime(2010, 1, 26)
__version__  = "0.1"
__status__   = "Development"

missingValue=-9.99e-35

stations1968=["results/IBM_1970_station_GeorgesBank.nc"]
stations1993=["results/IBM_1970_station_GeorgesBank.nc"]
stationNames=["Georges Bank"]

co=["blue","red"]
years=[1968,1993]
stName=0
 
#for station1968, station1993 in zip(stations1968,stations1993):
for station in stations1968:       
    fig=plt.figure(1)
    prey=0
        
    cdf=Dataset(station,"r")
    ax = fig.add_subplot(2,1,1)
     
    temp  =cdf.variables["temp"][:,:,:,:]
    
    time  =cdf.variables["time"][:]
    timeIndex=cdf.variables["timeIndex"][:,:,:]
    
    refDate=datetime.datetime(1948,1,1,0,0,0)
    time2=[]
    for t in range(len(time)):
        if time[t]<1000000:
            time2.append(refDate + datetime.timedelta(hours=time[t]))
    print "Start time %s and end time %s"%(time2[0],time2[-1])
    maxInd=len(time2)
    
    sgr      =cdf.variables["sgr"][:,:,:,:]
    depth   =cdf.variables["depth"][:,:,:,:]
    survival =cdf.variables["survival_probability"][:,:,:,:]
    dh       =cdf.variables["deltaH"][:]
     
    Ncohorts=len(depth[:,0,0,0])
    Nindividuals=len(depth[0,:,0,0])
    print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
    for cohort in range(Ncohorts):
        for individual in range(Nindividuals):
            ax.plot(time[0:maxInd],-depth[cohort,individual,0:maxInd,0],linewidth = 3)
            ax.set_ylim(-60, 0)
            title(str(stationNames[stName]))
                 
        ylabel('specific growth rate (%d$^{-1})$')
    gca().set_xticklabels([])
    
    ax = fig.add_subplot(2,1,2)
    for cohort in range(Ncohorts):
        for individual in range(Nindividuals):
            ax.plot(time[0:maxInd],sgr[cohort,individual,0:maxInd,0],linewidth = 3)
            #ax.set_ylim(-60, 0)
            #title(str(stationNames[stName]))
   # ax.set_title("station_"+str(os.path.basename(station[:-3]))+".png")
    gca().set_xticklabels([])
            
    timez=[]; datez=[]
    for dateInd in range(0,maxInd,1):
        t=refDate + datetime.timedelta(hours=time[dateInd])
        tt=str(t.month)+"/"+str(t.day)+"/"+str(t.hour) +"/"+str(t.min)
        datez.append(tt)
        timez.append(time[dateInd])
            
    xticks(timez, datez, rotation=-90)
    
    show()
  