from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 9, 9)
__modified__ = datetime.datetime(2010, 2, 18)
__version__  = "0.1"
__status__   = "Development, modified 9.9.2009, 18.02.2010"

missingValue=-9.99e-35

stationsCold=["results/IBM_2002_ESM_2001_2100_northsea.nc"]
stationsWarm=["results/IBM_2002_ESM_2001_2100_northsea.nc"]

stationNames=["North Sea"]

co=["blue","red"]
years=['Cold','Warm']
stName=0
stations=0
fig=plt.figure()
datez=[]
datezEmpty=[]

for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
    prey=0
           
    subplot(2,2,stations+1)
    for station in st:
        cdf=Dataset(station,"r")
         
        depth =cdf.variables["depth"][:,:,:]
        time  =cdf.variables["time"][:]
        timeIndex=cdf.variables["timeIndex"][:,:,:]
        
        refDate=datetime.datetime(2001,1,1,0,0,0)
        time2=[]
        for t in range(len(time)):
            if time[t]<1000000:
                time2.append(refDate + datetime.timedelta(hours=time[t]))
        print "Start time %s and end time %s"%(time2[0],time2[-1])
        maxInd=len(time2)
        
        sgr      =cdf.variables["sgr"][:,:,:]
        wgt      =cdf.variables["wgt"][:,:,:]
        deltaH   =cdf.variables["deltaH"][:]
        #print 'Manually set the number of cohorts: line 71'
        Ncohorts=len(depth[:,0,0]) -1
        Nindividuals=len(depth[0,:,0])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
        for cohort in range(Ncohorts):
            index1=int(timeIndex[cohort,0,0])+1
            index2=int(timeIndex[cohort,0,1])-1
            days=int((index2-index1)/(24.0/float(deltaH)))
         
            meanDEPTH24H=np.zeros((days),dtype=float64)
            rangeDEPTH24H=np.zeros((days,2),dtype=float64)
            meanIndDEPTH24H=np.zeros((Nindividuals,days),dtype=float64)
            stdDEPTH24H=np.zeros((days),dtype=float64)
            stdPlus=np.zeros((days),dtype=float64)
            stdMinus=np.zeros((days),dtype=float64)
          
            if cohort==0 and stNumber==0:
                time24H=np.zeros((Ncohorts,days),dtype=float64)
            
            minimum=99999
            maximum=-99999
            
            t=refDate + datetime.timedelta(hours=time[index1])
            tt=str(t.month)+"/"+str(t.day)
            datez.append(tt)
            datezEmpty.append("")
                
            for i in range(days):
                k1=(index1+i*24./float(deltaH)) -1
                k2=(index1+(i+1)*24./float(deltaH)) -1
                
                for j in range(Nindividuals):
                    meanIndDEPTH24H[j,i]=(mean(depth[cohort,j,k1:k2]))
                    for kkk in range(int(k2-k1)):
                        index=int(k1+kkk)
                        if depth[cohort,j,index] > maximum:
                            maximum = depth[cohort,j,index] 
                            rangeDEPTH24H[i,0]=maximum
                        if depth[cohort,j,index] < minimum:
                            minimum = depth[cohort,j,index] 
                            rangeDEPTH24H[i,0]=minimum
                            
                minimum=99999
                maximum=-99999
                if stNumber==0:
                    time24H[cohort,i]=time[k1]
               
                meanDEPTH24H[i]=np.mean(meanIndDEPTH24H[:,i])
                stdDEPTH24H[i]=std(meanIndDEPTH24H[:,i])
                if stNumber==0:
                    time24H[cohort,i]=time[k1]
            #
            plot(time24H[cohort,:],-meanDEPTH24H,color=co[stNumber],linewidth = 3, alpha=0.5)
            devNEG=np.zeros(meanDEPTH24H.shape)
            devPOS=np.zeros(meanDEPTH24H.shape)
            #
            for d in range(len(meanDEPTH24H)):
                devNEG[d] =min(0.0,-meanDEPTH24H[d]-stdDEPTH24H[d])
                devPOS[d] =min(0.0,-meanDEPTH24H[d]+stdDEPTH24H[d])
            fill_between(time24H[cohort,:], devPOS, devNEG, facecolor=co[stNumber],alpha=0.3)
            
            plot(time24H[cohort,0],-np.mean(meanDEPTH24H),color='black')
            plot(time24H[cohort,0],-np.mean(meanDEPTH24H),color=co[stNumber],marker='s')

            ylabel("Depth (m)")
            gca().set_xticklabels([])
            #title(str(stationNames[stName]))
         
            
        stNumber+=1
        annotate(str(stationNames[stations]), xy=(time24H[0,0],-40),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))   
    ax=gca()
    axis('tight')
    ax.set_ylim(-50, 0)
    ax.axhline(y=0.0, linewidth=2, color='grey')
    if stations in [0,1]:
        xticks(time24H[:,0], datezEmpty, rotation=-90)
    else:
        xticks(time24H[:,0], datez, rotation=-90)

    stations+=1

plotfile="results/All_stations_depth.png"
#plt.savefig(plotfile)
print "Saving to file: %s"%(plotfile)
stName+=1

show()