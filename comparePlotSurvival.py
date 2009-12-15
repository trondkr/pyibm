#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 9, 9)
__modified__ = datetime.datetime(2009, 9, 10)
__version__  = "0.1"
__status__   = "Development"

missingValue=-9.99e-35

stations1968=["results/IBM_1993_station_EastandWestGreenland.nc",
              "results/IBM_1993_station_Lofoten.nc",
              "results/IBM_1993_station_GeorgesBank.nc",
              "results/IBM_1993_station_NorthSea.nc",
              "results/IBM_1993_station_Iceland.nc"]

stations1993=["results/IBM_1993_station_EastandWestGreenland.nc",
              "results/IBM_1993_station_Lofoten.nc",
              "results/IBM_1993_station_GeorgesBank.nc",
              "results/IBM_1993_station_NorthSea.nc",
              "results/IBM_1993_station_Iceland.nc"]

stationNames=["East Greenland","Lofoten","Georges Bank","North Sea","Iceland"]

co=["blue","red"]
years=[1968,1993]
stName=0

for station1968, station1993 in zip(stations1968,stations1993):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(station1968)
    print "2->%s\n"%(station1993)
    st=[]
    st.append(station1968)
    st.append(station1993)
    stNumber=0
    
    fig=plt.figure()
    prey=0
        
    for station in st:
        cdf=Dataset(station,"r")
        ax = fig.add_subplot(2,1,1)
         
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
        wgt      =cdf.variables["wgt"][:,:,:,:]
        length   =cdf.variables["length"][:,:,:,:]
        survival =cdf.variables["survival_probability"][:,:,:,:]
    
        Ncohorts=len(depth[:,0,0,0])
        Nindividuals=len(depth[0,:,0,0])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
        for cohort in range(Ncohorts):
            index1=int(timeIndex[cohort,0,0])+1
            index2=int(timeIndex[cohort,0,1])-1
            days=int((index2-index1)/24)
            
            meanDEPTH24H=np.zeros((days),dtype=float64)
            rangeDEPTH24H=np.zeros((days,2),dtype=float64)
            meanIndDEPTH24H=np.zeros((Nindividuals,days),dtype=float64)
            stdDEPTH24H=np.zeros((days),dtype=float64)
            stdPlus=np.zeros((days),dtype=float64)
            stdMinus=np.zeros((days),dtype=float64)
            days=int((index2-index1)/24)
            if cohort==0 and stNumber==0:
                time24H=np.zeros((Ncohorts,days),dtype=float64)
            
            minimum=99999
            maximum=-99999
            
            for i in range(days):
                k1=index1+i*24
                k2=index1+(i+1)*24
               
                j=0
                for kk in range(24):
                    
                    for j in range(Nindividuals):
                        meanIndDEPTH24H[j,i]=(mean(depth[cohort,j,k1:k2,prey]))
                        
                        index=k1+kk
                        if depth[cohort,j,index,prey] > maximum:
                            maximum = depth[cohort,j,index,prey] 
                            rangeDEPTH24H[i,0]=maximum
                        if depth[cohort,j,index,prey] < minimum:
                            minimum = depth[cohort,j,index,prey] 
                            rangeDEPTH24H[i,0]=minimum
                            
                minimum=99999
                maximum=-99999
                if stNumber==0:
                    time24H[cohort,i]=time[k1]
                
                meanDEPTH24H[i]=np.mean(meanIndDEPTH24H[:,i])
                stdDEPTH24H[i]=std(meanIndDEPTH24H[:,i])
                
            ax.plot(time24H[cohort,:],-meanDEPTH24H,color=co[stNumber],linewidth = 3, alpha=0.5)
            ax.fill_between(time24H[cohort,:], -meanDEPTH24H+stdDEPTH24H, -meanDEPTH24H-stdDEPTH24H, facecolor=co[stNumber],alpha=0.3)
            ylabel("Depth (m)")
            gca().set_xticklabels([])
            title(str(stationNames[stName]))
             
        ax = fig.add_subplot(2,1,2)
     
        
        if stNumber==0:
            timePsur1=np.ones((Ncohorts),dtype=float64)
            timePsur2=np.ones((Ncohorts),dtype=float64)
        larvaPsur=np.ones((Ncohorts,Nindividuals),dtype=float64)
        meanPsur=np.ones((Ncohorts),dtype=float64)
        stdPsur=np.ones((Ncohorts),dtype=float64)
        fitness=np.ones((Ncohorts),dtype=float64)
        datez=[]
    
        for cohort in range(Ncohorts):
                
            index1=int(timeIndex[cohort,0,0])+1
            index2=int(timeIndex[cohort,0,1])-1
            widthIndex=index2-index1
            for ind in range(Nindividuals):
                larvaPsur[cohort,ind]=survival[cohort,ind,index2,prey]
                index=(index1+index2)/2
                if stNumber==0:
                    timePsur1[cohort]=time[index1]
                if stNumber==1:
                    timePsur2[cohort]=timePsur1[cohort]+widthIndex/2
                    
            t=refDate + datetime.timedelta(hours=time[index1])
            
            tt=str(t.month)+"/"+str(t.day)
            datez.append(tt)
           
            meanPsur[cohort]=np.mean(larvaPsur[cohort,:])
            stdPsur[cohort]=np.std(larvaPsur[cohort,:])
            widthIndex=index2-index1
            fitness[cohort]=meanPsur[cohort]*(mean(wgt[cohort,:,index2,prey],0))
            
            if cohort==Ncohorts-1 and stNumber==1:
                tline=[timePsur1[5] for hh in range(2)]       
                xline90=[0.04 for hh in range(2)]
                ax.plot(tline,xline90,color=co[1],linewidth = 4)
                ax.annotate("1993", xy=(timePsur1[6],0.14),  xycoords='data',
                    size=16,
                    bbox=dict(boxstyle="round", fc=co[stNumber], alpha=0.8))
                
               
            if cohort==Ncohorts-1 and stNumber==0:
                tline=[timePsur1[5] for hh in range(2)]       
                xline68=[0.03 for hh in range(2)] 
                ax.plot(tline,xline68,color=co[0],linewidth = 4)
                ax.annotate("1968", xy=(timePsur1[6],0.17),  xycoords='data',
                    size=16,
                    bbox=dict(boxstyle="round", fc=co[stNumber], alpha=0.5))
            
            
            print "Total survival %s for cohort is %s"%(meanPsur[cohort]*100.,cohort)
        if stNumber==0:
            ax.bar(timePsur1, meanPsur*100., width=widthIndex/3, color=co[stNumber], alpha=0.5, yerr=stdPsur)
        if stNumber==1:
             ax.bar(timePsur2,meanPsur*100., width=widthIndex/3, color=co[stNumber], alpha=0.8, yerr=stdPsur)
        ax.set_ylim(0, 0.2)
        xticks(timePsur2, datez, rotation=-90)
        ylabel("Survival probability (%)")
        stNumber+=1
   
    plotfile="results/compare_"+str(stationNames[stName])+"_survival.png"
    plt.savefig(plotfile)
    print "Saving to file: %s"%(plotfile)
    stName+=1
     
    show()
    #sys.exit()
      