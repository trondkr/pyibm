from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates

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

co=["slategrey","darkgrey"]
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
    preyList=[1,2]
        
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
        wgt      =cdf.variables["wgt"][:,:,:,:]
        length   =cdf.variables["length"][:,:,:,:]
        survival =cdf.variables["survival_probability"][:,:,:,:]
        deltaH   =cdf.variables["deltaH"][:]
        Ncohorts=len(depth[:,0,0,0]) 
        Nindividuals=len(depth[0,:,0,0])
        print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
       
        
        if stNumber==0:
            timePsur1=np.ones((Ncohorts),dtype=float64)
            timePsur2=np.ones((Ncohorts),dtype=float64)
            larvaPsur=np.ones((Ncohorts,Nindividuals,Nprey),dtype=float64)
            meanPsur=np.ones((Ncohorts,Nprey),dtype=float64)
            stdPsur=np.ones((Ncohorts,Nprey),dtype=float64)
            fitness=np.ones((Ncohorts,Nprey),dtype=float64)
            datez=[]
      
        for prey in range(len(preyList)):
           # print "Calculating average values for prey level %s"%(preyList[prey])
              
            for cohort in range(Ncohorts):
                 
                index1=int(timeIndex[cohort,0,0])+1
                index2=int(timeIndex[cohort,0,1])-1
                days=int((index2-index1)/(24.0/float(deltaH)))
                widthIndex=index2-index1
                    
                for ind in range(Nindividuals):
                    larvaPsur[cohort,ind,prey]=survival[cohort,ind,index2,preyList[prey]]
                   
                    if stNumber==0 and prey==0:
                        timePsur1[cohort]=time[index2]
                    if stNumber==0 and prey==0:
                        timePsur2[cohort]=timePsur1[cohort]+widthIndex/4.
               
                t=refDate + datetime.timedelta(hours=time[index2])
                tt=str(t.month)+"/"+str(t.day)
                datez.append(tt)
             
                meanPsur[cohort,prey]=np.mean(larvaPsur[cohort,:,prey])
              
                stdPsur[cohort,prey]=np.std(larvaPsur[cohort,:,prey])
                fitness[cohort,prey]=meanPsur[cohort,prey]*(mean(wgt[cohort,:,index2,prey]))
            """Here you choose what to plot (var):"""
            var=np.squeeze(meanPsur[:,prey])*100.
           # var=np.squeeze(fitness[:,prey])*100.
                
           # print "Total survival %s for cohort is %s at time %s for prey level %s\n"%(meanPsur[cohort,prey]*100.,cohort,t,prey)
            if stNumber==0 and prey==0:
                #plot(timePsur1, np.squeeze(meanPsur[:,prey])*100., color=co[stNumber],linewidth = 3,alpha=1)
                #plot(timePsur1, np.squeeze(meanPsur[:,prey])*100., color=co[stNumber],marker=styles[prey],alpha=1)
                bar(timePsur1, var, width=widthIndex/4.5, color=co[stNumber], alpha=0.8, yerr=stdPsur[:,prey]*100.)
                coldSurvival=meanPsur[cohort,prey]*100.
            if stNumber==1 and prey==0:
                #plot(timePsur1, np.squeeze(meanPsur[:,prey])*100., color=co[stNumber], linewidth = 3,alpha=1)
                #plot(timePsur1, np.squeeze(meanPsur[:,prey])*100., color=co[stNumber], marker=styles[prey],alpha=1)
                bar(timePsur2, var, width=widthIndex/4.5, color=co[stNumber], alpha=0.8, yerr=stdPsur[:,prey]*100.)
                warmSurvival=meanPsur[cohort,prey]*100.
          
            """Plot the higher prey density values as shaded colors in the background"""
            if stNumber==0 and prey==1:
                bar(timePsur1, var, width=widthIndex/4.5, color=co[stNumber], alpha=0.3, yerr=stdPsur[:,prey]*100.)
                coldSurvival=meanPsur[cohort,prey]*100.
            if stNumber==1 and prey==1:
                bar(timePsur2, var, width=widthIndex/4.5, color=co[stNumber], alpha=0.3, yerr=stdPsur[:,prey]*100.)
                warmSurvival=meanPsur[cohort,prey]*100.
                
            print "Average survival for year %s at station %s for prey %s is %s"%(station[12:16],stationNames[stName],prey,np.mean(meanPsur[:,prey])*100.)
        stNumber+=1
    
    
            
    ax=gca()
    axis('tight')  
    ax.set_ylim(0, 0.3)
    ax.annotate(str(stationNames[stName]), xy=(timePsur1[1],0.27),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))
    
    xticks(timePsur2, datez, rotation=-90)
    ylabel("Survival probability (%)")
    #ylabel("Fitness")
   # plotfile="results/compare_"+str(stationNames[stName])+"_fitness.png"
   # plt.savefig(plotfile)
    
    plotfile="results/compare_"+str(stationNames[stName])+"_survival.pdf"
    plt.savefig(plotfile)
    
    print "Saving to file: %s"%(plotfile)
    stName+=1
       
   # for t in range(len(timePsur2)):
   #     print "Survival in month %s : %3.3f warm: %3.3f cold: %3.3f"%(t,(warmSurvival[t]/(coldSurvival[t] + 0.00000001))*100., warmSurvival[t], coldSurvival[t])
    #show()
