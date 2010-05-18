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
__modified__ = datetime.datetime(2010, 3, 2)
__version__  = "0.1"
__status__   = "Development, modified 9.9.2009, 18.02.2010, 2.3.2010"

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
years=['Cold','Warm']
stName=0
stations=0

fig=plt.figure()
clf()
 
for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
    
    styles=['o','s','']
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
        
        sgr      =cdf.variables["sgr"][:,:,:,:]
        deltaH   =cdf.variables["deltaH"][:]
        Ncohorts=len(depth[:,0,0,0]) 
        Nindividuals=len(depth[0,:,0,0])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
        
        if stNumber==0:
            timeSGR1=np.ones((Ncohorts),dtype=float64)
            timeSGR2=np.ones((Ncohorts),dtype=float64)
            larvaSGR=np.ones((Ncohorts,Nindividuals,Nprey),dtype=float64)
            meanSGR=np.ones((Ncohorts,Nprey),dtype=float64)
            stdSGR=np.ones((Ncohorts,Nprey),dtype=float64)
            datez=[]
            datezEmpty=[]
        for prey in range(len(preyList)):
            for cohort in range(Ncohorts):
                   
                index1=int(timeIndex[cohort,0,0])+1
                index2=int(timeIndex[cohort,0,1])-1
                days=int((index2-index1)/(24.0/float(deltaH)))
                widthIndex=index2-index1
             
                for ind in range(Nindividuals):
                    larvaSGR[cohort,ind,prey]=sum(sgr[cohort,ind,index1:index2,preyList[prey]])/days
                   
                    if stNumber==0 and prey==0:
                        timeSGR1[cohort]=time[index1]
                    if stNumber==1 and prey==0:
                        timeSGR2[cohort]=timeSGR1[cohort]+widthIndex/4.
               
                t=refDate + datetime.timedelta(hours=time[index1])
                tt=str(t.month)+"/"+str(t.day)
                datez.append(tt)
                datezEmpty.append("")
                
                meanSGR[cohort,prey]=np.mean(larvaSGR[cohort,:,prey])
                stdSGR[cohort,prey]=np.std(larvaSGR[cohort,:,prey])
            
            pos=meanSGR[:,prey] + stdSGR[:,prey]
            neg=meanSGR[:,prey] - stdSGR[:,prey]
                
            if stNumber==0 and prey==0:
                fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3) 
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],linewidth = 4,alpha=1)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=1)
                #bar(timeSGR1, np.squeeze(meanSGR[:,prey]), width=widthIndex/4.5, color=co[stNumber], alpha=1.0, yerr=stdSGR[:,prey])
                coldSurvival=meanSGR[cohort,prey]*100.
            if stNumber==1 and prey==0:
                fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 4,alpha=1)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], marker=styles[prey],alpha=1)
                
                #bar(timeSGR2, np.squeeze(meanSGR[:,prey]), width=widthIndex/4.5, color=co[stNumber], alpha=1.0, yerr=stdSGR[:,prey])
                warmSurvival=meanSGR[cohort,prey]*100.
                
            """Plot the higher prey density values as shaded colors in the background"""
            if stNumber==0 and prey==1:
                fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 2,alpha=0.3)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=0.3)
                   
                #bar(timeSGR1, np.squeeze(meanSGR[:,prey]), width=widthIndex/4.5, color=co[stNumber], alpha=0.2, yerr=stdSGR[:,prey])
                coldSurvival=meanSGR[cohort,prey]*100.
            if stNumber==1 and prey==1:
                fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3) 
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 2,alpha=0.3)
                plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=0.3)
                #bar(timeSGR2, np.squeeze(meanSGR[:,prey]), width=widthIndex/4.5, color=co[stNumber], alpha=0.2, yerr=stdSGR[:,prey])
                warmSurvival=meanSGR[cohort,prey]*100.
        stNumber+=1
    ax=gca()
    axis('tight')
    ax.set_ylim(-4, 16)
    ax.axhline(y=0.0, linewidth=2, color='grey')
    if stations in [0,1]:
        xticks(timeSGR1, datezEmpty, rotation=-90)
    else:
        xticks(timeSGR1, datez, rotation=-90)
    ax.annotate(str(stationNames[stations]), xy=(timeSGR1[1],14),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))

    stations+=1    

#ylabel("Specific growth rate (%/day)")
plotfile="results/compare_"+str(stationNames[stName])+"_sgr.pdf"
plt.savefig(plotfile)
print "Saving to file: %s"%(plotfile)
stName+=1

show()
