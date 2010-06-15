from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

"""Import moduels to calculate daylength"""
sys.path.append('/Users/trond/Projects/arcwarm/light')
import calclight

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 9, 9)
__modified__ = datetime.datetime(2010, 6, 9)
__version__  = "0.1"
__status__   = "Development, modified 9.9.2009, 18.02.2010, 2.3.2010, 9.6.2010"

missingValue=-9.99e-35

stationsCold=["results/IBM_1977_station_Lofoten.nc",
              "results/IBM_1963_station_NorthSea.nc",
              "results/IBM_1983_station_Iceland.nc",
              "results/IBM_1965_station_GeorgesBank.nc"]

stationsWarm=["results/IBM_1990_station_Lofoten.nc",
              "results/IBM_1999_station_NorthSea.nc",
              "results/IBM_1960_station_Iceland.nc",
              "results/IBM_1999_station_GeorgesBank.nc"]

#stationsCold=["results/IBM_1977_station_Lofoten.nc"]

#stationsWarm=["results/IBM_1990_station_Lofoten.nc"]

latList=[67.5001, 54.5601, 63.7010,   41.6423]    
stationNames=["Lofoten","North Sea","Iceland","Georges Bank"]

co=["blue","red"]
co=["#0000cc","#FF0000"]
years=['Cold','Warm']
stName=0
stations=0

for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing station for two different years :"
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0
        
    fig=plt.figure()
    left, width = 0.1, 0.8
    rect1 = [left, 0.72, width, 0.20]
    rect2 = [left, 0.1, width, 0.60]
    ax1 = fig.add_axes(rect1)
    ax2 = fig.add_axes(rect2, sharex=ax1)

    styles=['o','s','']
    preyList=[0]
        
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
            datezEmpty=[]
            dayLength=np.ones((Ncohorts),dtype=np.float64)
            
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
               
                tt=t.timetuple()
                """Calculate daylength for given day"""
                radfl0=0.0; radmax=0.0; cawdir=0.0; clouds=0.0
                lat=latList[stations]
                radfl0,radmax,cawdir = calclight.calclight.qsw(radfl0,radmax,cawdir,clouds, lat*np.pi/180.0,int(tt[7])*1.0,365)
                radmax = radmax/0.217 # From W/m2 to umol/m2/s-1
                
                """Now calculate the total irradiance and daylength for given day"""
                daylen=0.0; totIrradiance = 0.
                diffAtt=0.18
                forageDepth=10.0
                dayLightThresh = 1
                dt_per_Day=960
                dh = 24./(dt_per_Day*1.0)	#Hours per timestep
                for hour in range(dt_per_Day):
                    h = dh*1.0*hour
                    sunHeight=0.0; surfLight=0.0
                    sunHeight, surfLight = calclight.calclight.surlig(h,radmax,int(tt[7])*1.0,lat,sunHeight,surfLight)
                    Eb = surfLight*exp(-diffAtt*forageDepth)
        
                    if Eb > dayLightThresh:
                     daylen = daylen + dh
                
                dayLength[cohort]=daylen
                meanSGR[cohort,prey]=np.mean(larvaSGR[cohort,:,prey])
                stdSGR[cohort,prey]=np.std(larvaSGR[cohort,:,prey])
                
            pos=meanSGR[:,prey] + stdSGR[:,prey]
            neg=meanSGR[:,prey] - stdSGR[:,prey]
                
            if stNumber==0 and prey==0:
                """Plot the growth rate at lower panel"""
                #ax2.fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color='k',linewidth = 5,alpha=1)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],linewidth = 4,alpha=1)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=1)
                
            if stNumber==1 and prey==0:
                #ax2.fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color='k',linewidth = 5,alpha=1)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 4,alpha=1)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], marker=styles[prey],alpha=1)
                
            """Plot the higher prey density values as shaded colors in the background"""
            if stNumber==0 and prey==1:
                #ax2.fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 2,alpha=0.3)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=0.3)
                
            if stNumber==1 and prey==1:
                #ax2.fill_between(timeSGR1,neg,pos, facecolor='grey',alpha=0.3) 
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber], linewidth = 2,alpha=0.3)
                ax2.plot(timeSGR1, np.squeeze(meanSGR[:,prey]), color=co[stNumber],marker=styles[prey],alpha=0.3)
            
        ax2.set_ylim(-4, 16)
        ax2.set_xlim(timeSGR1[0],timeSGR1[-1])
        ax1.set_xticklabels([])
     
        """Plot the daylength at the top panel"""
        ax1.fill_between(timeSGR1, 0, dayLength, color='#A5A5A5',linewidth = 0.1,alpha=1)     
       # ax1.plot(timeSGR1, dayLength, color='black',linewidth = 1,alpha=1)        
        ax1.set_ylim(0, 24.1)
        ax1.set_xlim(timeSGR1[0],timeSGR1[-1]) 
        ax1.set_xticklabels([""],color="grey", alpha=0.0)
    
        stNumber+=1
   
    ax2.axhline(y=0.0, linewidth=2, color='grey')
    xticks(timeSGR1, datez, rotation=-90)
    ax2.annotate(str(stationNames[stations]), xy=(timeSGR1[1],14),  xycoords='data',size=16, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))
    ax2.set_ylabel("Specific growth rate (%/day)")
    ax1.set_ylabel("Daylength (hours)")
    
    plotfile="results/sgr_"+str(stationNames[stations])+".pdf"
    plt.savefig(plotfile)
    print "Saving to file: %s"%(plotfile)
    
    stations+=1
    clf()
stName+=1