from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 9, 9)
__modified__ = datetime.datetime(2010, 4, 29)
__version__  = "0.1"
__status__   = "Development, modified 9.9.2009, 18.02.2010, 2.3.2010, 29.04.2010"

"""This script plots the nauplii concentration experienced by
a larva at a fixed depth through the year. Input files are created running ibm.py for cohorts of larva
that have been manually held at 10m depth."""

missingValue=-9.99e-35

stationsCold=["results/IBM_1977_station_Lofoten.nc",
              "results/IBM_1963_station_NorthSea.nc",
              "results/IBM_1983_station_Iceland.nc",
              "results/IBM_1965_station_GeorgesBank.nc"]

stationsWarm=["results/IBM_1990_station_Lofoten.nc",
              "results/IBM_1999_station_NorthSea.nc",
              "results/IBM_1960_station_Iceland.nc",
              "results/IBM_1999_station_GeorgesBank.nc"]

stationNames=["Lofoten","North Sea","Iceland","Georges Bank"]

co=["blue","red"]
colight=["blue","red"]
years=['Cold','Warm']
stName=0
stations=0
fig=plt.figure()
for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
    
    styles=['o','s','v','D']
    linestyles=[ '-' , '--' , '-.' , ':' ]
    preyList=[0]
    Nprey=len(preyList)
    subplot(2,2,stations+1)

    for station in st:
        cdf=Dataset(station,"r")
         
        depth =cdf.variables["depth"][:,:,:,:]
        time  =cdf.variables["time"][:]
        timeIndex=cdf.variables["timeIndex"][:,:,:]
        
        refDate=datetime.datetime(1948,1,1,0,0,0)
        time2=[]
        for t in range(len(time)):
            if time[t]<1000000:
                time2.append(refDate + datetime.timedelta(hours=time[t]))
                
        print "Start time %s and end time %s"%(time2[0],time2[-1])
        maxInd=len(time2)
        
        nauplii      =cdf.variables["nauplii"][:,:,:,:]
        deltaH   =cdf.variables["deltaH"][:]
        Ncohorts=len(depth[:,0,0,0]) 
        Nindividuals=len(depth[0,:,0,0])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
        
        if stNumber==0 and stations==0:
            timeNauplii1=np.ones((Ncohorts),dtype=float64)
            timeNauplii2=np.ones((Ncohorts),dtype=float64)
            larvaNauplii=np.ones((Ncohorts,Nindividuals,Nprey),dtype=float64)
            datez=[]; datezEmpty=[]
       
        for prey in range(len(preyList)):
            for cohort in range(Ncohorts):
                   
                index1=int(timeIndex[cohort,0,0])+1
                index2=int(timeIndex[cohort,0,1])-1
                days=int((index2-index1)/(24.0/float(deltaH)))
                widthIndex=index2-index1
             
                larvaNauplii[cohort,prey]=nauplii[cohort,0,index1,preyList[prey]]
                
                if stNumber==0 and prey==0 and stations==0:
                    timeNauplii1[cohort]=time[index1]
                    t=refDate + datetime.timedelta(hours=time[index1])
                    tt=str(t.month)+"/"+str(t.day)
                    datez.append(tt)
                    datezEmpty.append("")
                
                if stNumber==1 and prey==0 and stations==0:
                    timeNauplii2[cohort]=timeNauplii1[cohort]+widthIndex/4.
               
                
              
            if stNumber==0 and prey==0:
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=co[stNumber],linewidth = 3,alpha=1)
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=co[stNumber],marker=styles[0],alpha=1)
            if stNumber==0 and prey==1:
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=colight[stNumber],linewidth = 3,alpha=0.5)
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=colight[stNumber],marker=styles[1],alpha=0.5)
                
            if stNumber==1 and prey==0:
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=co[stNumber], linewidth = 3,alpha=1)
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=co[stNumber], marker=styles[2],alpha=1)
            if stNumber==1 and prey==1:
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=colight[stNumber], linewidth = 3,alpha=0.5)
                plot(timeNauplii1, np.squeeze(larvaNauplii[:,prey]), color=colight[stNumber], marker=styles[3],alpha=0.5)
            
        stNumber+=1
    ax=gca()
    axis('tight')
    ax.set_ylim(0, 30)
    ax.axhline(y=0.0, linewidth=2, color='grey')
    if stations in [0,1]:
        xticks(timeNauplii1, datezEmpty, rotation=-90)
    else:
        xticks(timeNauplii1, datez, rotation=-90)
    ax.annotate(str(stationNames[stations]), xy=(timeNauplii1[1],26),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))
  
    print "Plotting station %s"%(stationNames[stations])
    stations+=1
    
 

#ylabel("Specific growth rate (%/day)")
plotfile="results/All_stations_Nauplii.pdf"
plt.savefig(plotfile)
print "Saving to file: %s"%(plotfile)
stName+=1

show()
