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
              "results/IBM_1965_station_GeorgesBank.nc",
              "results/IBM_1963_station_NorthSea.nc",
              "results/IBM_1983_station_Iceland.nc"]

stationsWarm=["results/IBM_1990_station_Lofoten.nc",
              "results/IBM_1999_station_GeorgesBank.nc",
              "results/IBM_1999_station_NorthSea.nc",
              "results/IBM_1960_station_Iceland.nc"]

stationNames=["Lofoten","Georges Bank","North Sea","Iceland"]

co=["blue","red"]
years=['Cold','Warm']
stName=0

fig=plt.figure()
for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
    clf()    
    
    styles=['o','s','']
    preyList=[0,1]
        
    Nprey=len(preyList)
    
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
    ax.set_ylim(-4, 20)
    ax.axhline(y=0.0, linewidth=2, color='grey')        
    ax.annotate(str(stationNames[stName]), xy=(timeSGR1[1],17),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))
  
    xticks(timeSGR2, datez, rotation=-90)
    ylabel("Specific growth rate (%/day)")
    plotfile="results/compare_"+str(stationNames[stName])+"_sgr.png"
    plt.savefig(plotfile)
    print "Saving to file: %s"%(plotfile)
    stName+=1
       
    # for t in range(len(timeSGR2)):
    #     print "Survival in month %s : %3.3f warm: %3.3f cold: %3.3f"%(t,(warmSurvival[t]/(coldSurvival[t] + 0.00000001))*100., warmSurvival[t], coldSurvival[t])
    show()
