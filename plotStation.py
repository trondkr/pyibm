#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

years    = mdates.YearLocator()   # every year
months   = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 8, 13)
__modified__ = datetime.datetime(2009, 8, 13)
__version__  = "0.1"
__status__   = "Development"

missingValue=-9.99e-35

stations=["results/IBM_1971_station_EastandWestGreenland.nc",
         "results/IBM_1971_station_Lofoten.nc",
         "results/IBM_1971_station_GeorgesBank.nc",
         "results/IBM_1971_station_NorthSea.nc",
         "results/IBM_1971_station_Iceland.nc"]

for station in stations:
    cdf=Dataset(station,"r")
    
    depth=cdf.variables["depth"][:,:,:,:]
    time=(cdf.variables["time"][:])
    timeIndex=(cdf.variables["timeIndex"][:,:,:])
    
    refDate=datetime.datetime(1948,1,1,0,0,0)
    time2=[]
    for t in range(len(time)):
        if time[t]<1000000:
            time2.append(refDate + datetime.timedelta(hours=time[t]))
    print "Start time %s and end time %s"%(time2[0],time2[-1])
    maxInd=len(time2)
    
    sgr=cdf.variables["sgr"][:,:,:,:]
    wgt=cdf.variables["wgt"][:,:,:,:]
    length=cdf.variables["length"][:,:,:,:]
    survival=cdf.variables["survival_probability"][:,:,:,:]
    
    sgrrel=cdf.variables["sgr_rel"][:,:,:,:]
    
    
    fig=plt.figure()
    ax = fig.add_subplot(4,1,1)
    colors=('r','g','b','c','y','m','k','r','g','b','c','y','m','k','r','g','b','c','y','m','k')
    
    """Calculate 24 hour running mean of specific growth rate"""
    Ker = np.ones((48,), dtype='float')
    Ker[0], Ker[-1] = 0.5, 0.5
    Ker = Ker/(48.)
    depth24H=np.zeros(depth.shape,dtype=float64)
    
    prey=0
    Ncohorts=len(depth[:,0,0,0])
    Nindividuals=len(depth[0,:,0,0])
    print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
    for cohort in range(Ncohorts):
        index1=int(timeIndex[cohort,0,0])+1
        index2=int(timeIndex[cohort,0,1])-1
        if index2==missingValue:
            index2=len(timeIndex-1)
       
        days=int((index2-index1)/24)
        #print "Finding %s days period of growths for cohort %s"%(days,cohort)
        
        meanSGR24H=np.zeros((days),dtype=float64)
        rangeSGR24H=np.zeros((days,2),dtype=float64)
        meanIndSGR24H=np.zeros((Nindividuals,days),dtype=float64)
        stdSGR24H=np.zeros((days),dtype=float64)
        stdPlus=np.zeros((days),dtype=float64)
        stdMinus=np.zeros((days),dtype=float64)
        
        time24H=np.zeros((days),dtype=float64)
        
        minimum=99999
        maximum=-99999
        
        for i in range(days):
            k1=index1+i*24
            k2=index1+(i+1)*24
           
            meanSGR24H[i]=sum(sgr[cohort,:,k1:k2,prey])/Nindividuals
            j=0
            for j in range(Nindividuals):
                meanIndSGR24H[j,i]=sum(sgr[cohort,j,k1:k2,prey])
                
                if sum(sgr[cohort,j,k1:k2,prey])> maximum:
                    maximum =  sum(sgr[cohort,j,k1:k2,prey])
                if sum(sgr[cohort,j,k1:k2,prey])< minimum:
                    minimum =  sum(sgr[cohort,j,k1:k2,prey])
                    
            rangeSGR24H[i,0]=minimum
            rangeSGR24H[i,1]=maximum
            minimum=99999
            maximum=-99999
        
            stdSGR24H[i]=std(meanIndSGR24H[:,i])
            time24H[i]=time[k1]
            
        for i in range(days):
            stdMinus[i] = meanSGR24H[i] - stdSGR24H[i]
            stdPlus[i]  = meanSGR24H[i] + stdSGR24H[i]
            
     #   sgr24H[cohort,ind,index1:index2,0]= np.convolve(sgr[cohort,ind,index1:index2,0], Ker, mode='same')
        ax.plot(time24H,meanSGR24H,color="black",linewidth = 3)
        ax.fill_between(time24H, stdMinus, stdPlus, facecolor='black',alpha=0.4)
        ax.fill_between(time24H, rangeSGR24H[:,0], rangeSGR24H[:,1], facecolor='black',alpha=0.2)
        
        #ax.plot(time[index1:index2],sgr[cohort,ind,index1:index2,0],colors[cohort],linewidth = 2)
        
    ax.set_title("station_"+str(os.path.basename(station[:-3]))+".png")
    gca().set_xticklabels([])
    
    ax = fig.add_subplot(4,1,2)
     
    for cohort in range(Ncohorts):
    
        index1=int(timeIndex[cohort,0,0])+1
        index2=int(timeIndex[cohort,0,1])-1
        meanWGT = mean(wgt[cohort,:,index1:index2,prey],0)
        rangeWGT=np.zeros((len(meanWGT),2),dtype=float64)
        minimum=99999
        maximum=-99999
        i=0
        for i in range(index2-index1):
            j=i+index1
            if np.max(wgt[cohort,:,j,prey]) > maximum:
                maximum =  np.max(wgt[cohort,:,j,prey])
                rangeWGT[i,1]=maximum
            if np.min(wgt[cohort,:,j,prey]) < minimum:
                minimum = np.min(wgt[cohort,:,j,prey])
                rangeWGT[i,0]=minimum    
        
            minimum=99999
            maximum=-99999
       # print "Maximum %s and minimum %s and mean %s wgt and cohort %s"%(rangeWGT[i,1], rangeWGT[i,0],meanWGT.mean(), cohort)
       
        ax.plot(time[index1:index2],meanWGT,color="black",linewidth = 2)
        stdMinus=mean(wgt[cohort,:,index1:index2,prey],0) - std(wgt[cohort,:,index1:index2,prey],0)
        stdPlus=mean(wgt[cohort,:,index1:index2,prey],0) + std(wgt[cohort,:,index1:index2,prey],0)
        ax.fill_between(time[index1:index2], stdMinus, stdPlus, facecolor='black',alpha=0.4)
        ax.fill_between(time[index1:index2], rangeWGT[:,0], rangeWGT[:,1], facecolor='green',alpha=0.2)
        gca().set_xticklabels([])
        
    #ax.set_title('Weight')
    
    ax = fig.add_subplot(4,1,3)
    
    for cohort in range(Ncohorts):#len(depth[:,0,0,0])-1)
            
        index1=int(timeIndex[cohort,0,0])+1
        index2=int(timeIndex[cohort,0,1])-1
        meanDEPTH24H=np.zeros((days),dtype=float64)
        rangeDEPTH24H=np.zeros((days,2),dtype=float64)
        meanIndDEPTH24H=np.zeros((Nindividuals,days),dtype=float64)
        stdDEPTH24H=np.zeros((days),dtype=float64)
        stdPlus=np.zeros((days),dtype=float64)
        stdMinus=np.zeros((days),dtype=float64)
        
        time24H=np.zeros((days),dtype=float64)
        
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
        
            time24H[i]=time[k1]
            
            meanDEPTH24H[i]=np.mean(meanIndDEPTH24H[:,i])
            stdDEPTH24H[i]=std(meanIndDEPTH24H[:,i])
            
        ax.plot(time24H,-meanDEPTH24H,color="black",linewidth = 3)
        ax.fill_between(time24H, -meanDEPTH24H+stdDEPTH24H, -meanDEPTH24H-stdDEPTH24H, facecolor='green',alpha=0.4)
        
        gca().set_xticklabels([])
      #  ax.fill_between(time24H, -rangeDEPTH24H[:,1], -rangeDEPTH24H[:,0], facecolor='green',alpha=0.5)
      
    #ax.set_title('Depth')
    
    ax = fig.add_subplot(4,1,4)
    
    timePsur=np.ones((Ncohorts),dtype=float64)
    larvaPsur=np.ones((Ncohorts,Nindividuals),dtype=float64)
    meanPsur=np.ones((Ncohorts),dtype=float64)
    stdPsur=np.ones((Ncohorts),dtype=float64)
    fitness=np.ones((Ncohorts),dtype=float64)
    datez=[]

    for cohort in range(Ncohorts):#len(depth[:,0,0,0])-1)
            
        index1=int(timeIndex[cohort,0,0])+1
        index2=int(timeIndex[cohort,0,1])-1
        for ind in range(Nindividuals):
            larvaPsur[cohort,ind]=survival[cohort,ind,index2,prey]
            index=(index1+index2)/2
            timePsur[cohort]=time[index1]
       
        t=refDate + datetime.timedelta(hours=time[index1])
        
        tt=str(t.year)+"/"+str(t.month)+"/"+str(t.day)
        datez.append(tt)
       
        meanPsur[cohort]=np.mean(larvaPsur[cohort,:])
        stdPsur[cohort]=np.std(larvaPsur[cohort,:])
        widthIndex=index2-index1
        fitness[cohort]=meanPsur[cohort]*(mean(wgt[cohort,:,index2,prey],0))
      
        print "Total survival %s for cohort is %s"%(meanPsur[cohort]*100.,cohort)
    
    ax.bar(timePsur, meanPsur, width=widthIndex/2, color='green', alpha=0.5, yerr=stdPsur)
    
    xticks(timePsur, datez, rotation=-90)
#    ax.set_title('Depth')
    
    plotfile="station_"+str(os.path.basename(station[:-3]))+".png"
    plt.savefig(plotfile)
    print "Saving to file: %s"%(plotfile)
    
    #show()