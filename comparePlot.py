
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 8, 13)
__modified__ = datetime.datetime(2009, 8, 13)
__version__  = "0.1"
__status__   = "Development"

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
 
for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "Comparing station at two different years:"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
    
    fig=plt.figure()
    prey=0
        
    for station in st:
        cdf=Dataset(station,"r")
        ax = fig.add_subplot(2,1,1)
         
        depth =cdf.variables["depth"][:,:,:,:]
        temp  =cdf.variables["temp"][:,:,:,:]
        
        time  =cdf.variables["time"][:]
        timeIndex=cdf.variables["timeIndex"][:,:,:]
        
        refDate=datetime.datetime(1948,1,1,0,0,0)
        time2=[]
        for t in range(len(time)):
            if time[t]<1000000:
                time2.append(refDate + datetime.timedelta(hours=time[t]))
#                print "hour",datetime.timedelta(hours=time[t])
            
        maxInd=len(time2)
        
        sgr      =cdf.variables["sgr"][:,:,:,:]
        sgrAF    =cdf.variables["sgrAF"][:,:,:,:]
        length   =cdf.variables["length"][:,:,:,:]
        survival =cdf.variables["survival_probability"][:,:,:,:]
        deltaH   =cdf.variables["deltaH"][:]
         
        Ncohorts=len(depth[:,0,0,prey]) - 1
        Nindividuals=len(depth[0,:,0,prey])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
        for cohort in range(Ncohorts):
            index1=int(timeIndex[cohort,0,0])+1
            index2=int(timeIndex[cohort,0,1])-1
            if index2==missingValue:
                index2=len(timeIndex-1)
           
            days=int((index2-index1)/(24.0/deltaH))
        
            meanSGR24H      =np.zeros((days),dtype=float64)
            meanSGRAF24H      =np.zeros((days),dtype=float64)
            meanIndSGR24H   =np.zeros((Nindividuals,days),dtype=float64)
            meanIndSGRAF24H =np.zeros((Nindividuals,days),dtype=float64)
            stdSGR24H       =np.zeros((days),dtype=float64)
            stdSGRAF24H       =np.zeros((days),dtype=float64)
            stdPlus       =np.zeros((days),dtype=float64)
            stdMinus      =np.zeros((days),dtype=float64)
            rangeSGR24H   =np.zeros((days,2),dtype=float64)
            
            if stNumber == 0 and cohort==0:
                """Statistics:"""
                meanSGRWarm       =np.zeros((Ncohorts),dtype=float64)
                meanSGRAFWarm     =np.zeros((Ncohorts),dtype=float64)
                stdSGRAFWarm      =np.zeros((Ncohorts),dtype=float64)
                stdSGRWarm        =np.zeros((Ncohorts),dtype=float64)
                meanSGRCold       =np.zeros((Ncohorts),dtype=float64)
                meanSGRAFCold     =np.zeros((Ncohorts),dtype=float64)
                stdSGRAFCold      =np.zeros((Ncohorts),dtype=float64)
                stdSGRCold        =np.zeros((Ncohorts),dtype=float64)
                
                time24H       =np.zeros((Ncohorts,days),dtype=float64)
            
            minimum=99999
            maximum=-99999
            
            for i in range(days):
                k1=(index1+i*24./float(deltaH)) -1
                k2=(index1+(i+1)*24./float(deltaH)) -1
              
                meanSGR24H[i] = sum(sgr[cohort,:,k1:k2,prey])/Nindividuals
                meanSGRAF24H[i] = sum(sgrAF[cohort,:,k1:k2,prey])/Nindividuals
                
                j=0
                for j in range(Nindividuals):
    
                    meanIndSGR24H[j,i] = sum(sgr[cohort,j,k1:k2,prey])
                    meanIndSGRAF24H[j,i] = sum(sgrAF[cohort,j,k1:k2,prey])
                    
                    if sum(sgr[cohort,j,k1:k2,prey])> maximum:
                        maximum =  sum(sgr[cohort,j,k1:k2,prey])
                    if sum(sgr[cohort,j,k1:k2,prey])< minimum:
                        minimum =  sum(sgr[cohort,j,k1:k2,prey])
                        
                rangeSGR24H[i,0]=minimum
                rangeSGR24H[i,1]=maximum
                minimum=99999
                maximum=-99999
            
                stdSGR24H[i]   = std(meanIndSGR24H[:,i])
                stdSGRAF24H[i] = std(meanIndSGRAF24H[:,i])
            
                if stNumber == 0:    
                    time24H[cohort,i]   = time[k1]
                
            if stNumber==0:
                meanSGRCold[cohort]=np.mean(meanSGR24H)
                stdSGRCold[cohort]=np.mean(stdSGR24H)
                meanSGRAFCold[cohort]=np.mean(meanSGRAF24H)
                stdSGRAFCold[cohort]=np.mean(stdSGRAF24H)
              
            if stNumber==1:
                meanSGRWarm[cohort]=np.mean(meanSGR24H)
                stdSGRWarm[cohort]=np.mean(stdSGR24H)
                meanSGRAFWarm[cohort]=np.mean(meanSGRAF24H)
                stdSGRAFWarm[cohort]=np.mean(stdSGRAF24H)
               
            ax.plot(time24H[cohort,:],meanSGR24H,color=co[stNumber],linewidth = 3)
            ax.fill_between(time24H[cohort,:], stdMinus, stdPlus, facecolor=co[stNumber],alpha=0.4)
            ax.fill_between(time24H[cohort,:], rangeSGR24H[:,0], rangeSGR24H[:,1], facecolor=co[stNumber],alpha=0.2)
            ax.set_ylim(2, 25)
            title(str(stationNames[stName]))
              
            ylabel('specific growth rate (%d$^{-1})$')
        gca().set_xticklabels([])
       
      
        ax = fig.add_subplot(2,1,2)
        aveTCold=0.0; aveTWarm=0.0
        for cohort in range(Ncohorts):
        
            index1=int(timeIndex[cohort,0,0])+1
            index2=int(timeIndex[cohort,0,1])-1
           
            meanLENGTH = mean(length[cohort,:,index1:index2,prey],0)
            rangeLENGTH=np.zeros((len(meanLENGTH),2),dtype=float64)
            minimum=99999
            maximum=-99999
            i=0
            for i in range(index2-index1):
                j=i+index1
                if np.max(length[cohort,:,j,prey]) > maximum:
                    maximum =  np.max(length[cohort,:,j,prey])
                    rangeLENGTH[i,1]=maximum
                if np.min(length[cohort,:,j,prey]) < minimum:
                    minimum = np.min(length[cohort,:,j,prey])
                    rangeLENGTH[i,0]=minimum    
            
                minimum=99999
                maximum=-99999
           
            if cohort==0 and stNumber==0:
                saveIndex1=index1    
                oldtime=time
                datez=[]
                timez=[]
            if stNumber==0:
                t=refDate + datetime.timedelta(hours=time[index1])
                tt=str(t.month)+"/"+str(t.day)
                datez.append(tt)
                timez.append(time[index1])
                
            if stNumber==0:
                aveTCold=aveTCold + np.mean(np.mean(temp[cohort,:,index1:index2,prey],1))
                print "COLD: average T %s for cohort %s"%(np.mean(np.mean(temp[cohort,:,index1:index2,prey],1)),cohort)
                tcold="COLD (aveT=%2.2f)"%(aveTCold/Ncohorts)
            
            if stNumber==1:
                aveTWarm=aveTWarm + np.mean(np.mean(temp[cohort,:,index1:index2,prey],1))
                print "WARM: average T %s for cohort %s"%(np.mean(np.mean(temp[cohort,:,index1:index2,prey],1)),cohort)
                twarm="WARM (aveT=%2.2f)"%(aveTWarm/Ncohorts)
            
            #if cohort==Ncohorts-1 and stNumber==1:
            #    tline=[time24H[hh,0] for hh in range(2)]       
            #    xline90=[20.5 for hh in range(2)]
            #    ax.plot(tline,xline90,color=co[1],linewidth = 4)
            #    ax.annotate(t90, xy=(time24H[1,0],20),  xycoords='data',
            #        size=16,
            #        bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), alpha=1.0))
            #    
               
            #if cohort==Ncohorts-1 and stNumber==0:
            #    tline=[time24H[hh,0] for hh in range(2)]       
            #    xline68=[23.5 for hh in range(2)] 
            #    ax.plot(tline,xline68,color=co[0],linewidth = 4)
            #    ax.annotate(t68, xy=(time24H[1,0],23),  xycoords='data',
            #        size=16,
            #        bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), alpha=1.0))
            #
            ax.plot(oldtime[index1:index2],meanLENGTH,color=co[stNumber],linewidth = 2)
            stdMinus=mean(length[cohort,:,index1:index2,prey],0) - std(length[cohort,:,index1:index2,prey],0)
            stdPlus=mean(length[cohort,:,index1:index2,prey],0) + std(length[cohort,:,index1:index2,prey],0)
            ax.fill_between(oldtime[index1:index2], stdMinus, stdPlus, facecolor=co[stNumber],alpha=0.4)
            ax.fill_between(oldtime[index1:index2], rangeLENGTH[:,0], rangeLENGTH[:,1], facecolor=co[stNumber],alpha=0.2)
            ax.set_ylim(5, 25)
            ylabel("Length (mm)")
            xticks(timez, datez, rotation=-90)
           
            
        stNumber+=1   
   
    plotfile="results/compare_"+str(stationNames[stName])+"_length.png"
    plt.savefig(plotfile)
    print "Saving to file: %s\n"%(plotfile)
    stName+=1
     
    for cohort in range(len(datez)):
        print "Mean SGR in month %s : warm/cold %3.3f cold: %3.3f warm: %3.3f"%(datez[cohort],
                                                                                     100.*(np.mean(meanSGRWarm[cohort])/(np.mean(meanSGRCold[cohort])+0.00000001)),
                                                                                     np.mean(meanSGRCold[cohort]), np.mean(meanSGRWarm[cohort]))
        
    
    show()
   