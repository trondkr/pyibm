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
              "results/IBM_1963_station_NorthSea.nc",
              "results/IBM_1983_station_Iceland.nc",
              "results/IBM_1965_station_GeorgesBank.nc"]

stationsWarm=["results/IBM_1990_station_Lofoten.nc",
              "results/IBM_1999_station_NorthSea.nc",
              "results/IBM_1960_station_Iceland.nc",
              "results/IBM_1999_station_GeorgesBank.nc"]

stationNames=["Lofoten","North Sea","Iceland","Georges Bank"]

spawningStart=[3,1,4,2]
spawningEnd  =[5,5,5,5]
co=["slategrey","darkgrey"]
co=["blue","red"]

years=['Cold','Warm']
stName=0

"""Choose which variable to plot"""
variable="fitness"
variable="survival"

for stationCold, stationWarm in zip(stationsCold,stationsWarm):
    print "\n\nComparing %s at station for two different years :"%(variable)
    print "1->%s"%(stationCold)
    print "2->%s\n"%(stationWarm)
    st=[]
    st.append(stationCold)
    st.append(stationWarm)
    stNumber=0

    styles=['o','s','']
    preyList=[0]

    Nprey=len(preyList)
    fig=plt.figure()
    left, width = 0.1, 0.8
    rect1 = [left, 0.7, width, 0.25]
    rect2 = [left, 0.1, width, 0.60]
    ax1 = fig.add_axes(rect1)
    ax2 = fig.add_axes(rect2, sharex=ax1)
#    clf()
    for station in st:
        cdf=Dataset(station,"r")

        depth =cdf.variables["depth"][:,:,:,:]
        time  =cdf.variables["time"][:]
        timeIndex=cdf.variables["timeIndex"][:,:,:]
        length =cdf.variables["length"][:]

        print "Average length of individuals at all cohorts %3.4f+-%4.4f"%(np.mean(length),np.std(length))

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
            pdfCold=np.ones((Ncohorts,Nprey),dtype=float64)
            pdfWarm=np.ones((Ncohorts,Nprey),dtype=float64)
            stdPsur=np.ones((Ncohorts,Nprey),dtype=float64)
            fitness=np.ones((Ncohorts,Nprey),dtype=float64)
            spawning=np.ones((2),dtype=float64)
            datez=[]

        for prey in range(len(preyList)):
            for cohort in range(Ncohorts):

                index1=int(timeIndex[cohort,0,0])
                index2=int(timeIndex[cohort,0,1])-1
                days=int((index2-index1)/(24.0/float(deltaH)))
                widthIndex=((index2-index1))/4.1

                for ind in range(Nindividuals):
                    larvaPsur[cohort,ind,prey]=survival[cohort,ind,:,preyList[prey]].min()

                    if stNumber==0 and prey==0:
                        timePsur1[cohort]=time[index2]
                        if ind==0:
                            t=refDate +datetime.timedelta(hours=time[index2])

                            if spawningStart[stName]==int(t.month):
                                spawning[0]=cohort
                            if spawningEnd[stName]==int(t.month):
                                spawning[1]=cohort + 1
                    if ind==0:
                        print "Time : %s to %s (%s days) survival: %s"%(time2[index1],time2[index2], days,survival[cohort,ind,:,preyList[prey]].min())
                    if stNumber==0 and prey==0:
                        timePsur2[cohort]=timePsur1[cohort] +widthIndex


                t=refDate + datetime.timedelta(hours=time[index2])
                tt=str(t.day)+"/"+str(t.month)
                datez.append(tt)

                meanPsur[cohort,prey]=np.mean(larvaPsur[cohort,:,prey])

                stdPsur[cohort,prey]=np.std(larvaPsur[cohort,:,prey])
                fitness[cohort,prey]=meanPsur[cohort,prey]*(mean(wgt[cohort,:,-1,prey]))
            """Here you choose what to plot (var):"""
            if variable=="survival":
                var=np.squeeze(meanPsur[:,prey])*100.
            if variable=="fitness":
                var=np.squeeze(fitness[:,prey])*100.

            """Calculate the cumulative Probability Density Function and plot it"""
            if stNumber==0 and prey==0:
                for p in range(len(meanPsur[:,0])):
                    if variable=="fitness":
                        pdfCold[p,prey] = (np.sum(fitness[0:p,prey]))
                        totalCold=np.sum(np.sum(fitness[:,prey]))
                    if variable=="survival":
                        pdfCold[p,prey] = (np.sum(meanPsur[0:p,prey]))
                        totalCold=np.sum(np.sum(meanPsur[:,prey]))
            if stNumber==1 and prey==0:
                for p in range(len(meanPsur[:,0])):
                    if variable=="fitness":
                        pdfWarm[p,prey] = (np.sum(fitness[0:p,prey]))
                        totalWarm=np.sum(np.sum(fitness[:,prey]))
                    if variable=="survival":
                        pdfWarm[p,prey] = (np.sum(meanPsur[0:p,prey]))
                        totalWarm=np.sum(np.sum(meanPsur[:,prey]))

            if stNumber==0 and prey==0:
                ax2.bar(timePsur1, var, width=widthIndex, color=co[stNumber], alpha=0.8, yerr=stdPsur[:,prey]*100.)

            if stNumber==1 and prey==0:
                ax2.bar(timePsur2, var, width=widthIndex, color=co[stNumber], alpha=0.8, yerr=stdPsur[:,prey]*100.)

            """Plot the higher prey density values as shaded colors in the background"""
            if stNumber==0 and prey==1:
                ax2.bar(timePsur1, var, width=widthIndex, color=co[stNumber], alpha=0.3, yerr=stdPsur[:,prey]*100.)
            if stNumber==1 and prey==1:
                ax2.bar(timePsur2, var, width=widthIndex, color=co[stNumber], alpha=0.3, yerr=stdPsur[:,prey]*100.)


            ax2.set_ylim(0, 0.3)
            ax2.set_xlim(timePsur1[0],timePsur1[-1])
            ax1.set_xticklabels([])
            print "spawning start %s and end %s"%(datez[int(spawning[0])], datez[int(spawning[1])])
            if stName!=3:
                cumulative=(np.sum(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))*100.)
                print "================================================================================================="
                print "Cumulative %s during spawning for year %s at station %s for prey %s is %s"%(variable,
                                                                                                   station[12:16],
                                                                                                   stationNames[stName],
                                                                                                   prey,cumulative)
                print "Average survival rate was %3.4f +- %3.8f (standard error)"%(np.mean(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))*100.,
                                                                  np.std(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))/
                                                                  np.sqrt(len(meanPsur[spawning[0]:spawning[1],prey]))*100.)
                print "================================================================================================="
            else:
                #GB special case
                cumulative=(np.sum(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))*100.) + (np.sum(np.squeeze(meanPsur[21:-1,prey]))*100.)
                print "================================================================================================="
                print "Cumulative %s during spawning for year %s at station %s for prey %s is %s"%(variable,
                                                                                                   station[12:16],
                                                                                                   stationNames[stName],
                                                                                                   prey, cumulative)
              
                print "Average survival rate was %3.6f +- %3.8f (standard error)"%(np.mean(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))*100.+
                                                                  + (np.mean(np.squeeze(meanPsur[20:-1,prey]))*100.),
                                                                  (np.std(np.squeeze(meanPsur[spawning[0]:spawning[1],prey]))/
                                                                  np.sqrt(len(meanPsur[spawning[0]:spawning[1],prey]))*100.)+
                                                                  + (np.std(np.squeeze(meanPsur[20:-1,prey]))/
                                                                     np.sqrt(len(meanPsur[20:-1,prey]))*100.))
                
                print "================================================================================================="
        stNumber+=1

    ax1.plot(timePsur1, pdfWarm/(totalWarm+totalCold), color="red",linewidth = 3,alpha=1)
    ax1.plot(timePsur1, pdfCold/(totalWarm+totalCold), color="blue", linewidth = 3,alpha=1)

    ax1.set_ylim(0, 1)
    ax1.set_xlim(timePsur1[0],timePsur1[-1])
    ax1.set_xticklabels([""],color="grey", alpha=0.0)

    """Calculate the difference between a cold and a warm year in probability of survival"""
    print "Probability of survival in warm versus cold year is %s\n"%((np.mean(pdfWarm[spawning[0]:spawning[1],0])
                                                                    /np.mean(pdfCold[spawning[0]:spawning[1],0]))*100.)

    ax1.annotate(str(stationNames[stName]), xy=(timePsur1[1],0.8),  xycoords='data',size=16) #, bbox=dict(boxstyle="round", fc='grey', alpha=0.5))
    xticks(timePsur1,datez, rotation=-90)

    if variable=="survival":
        plotfile="results/compare_"+str(stationNames[stName])+"_survival.pdf"
        ax1.set_ylabel("PDF")
        ax2.set_ylabel("Survival prob. (%)")

    if variable=="fitness":
        ax2.set_ylabel("Fitness")
        ax1.set_ylabel("PDF")
        plotfile="results/compare_"+str(stationNames[stName])+"_fitness.pdf"

    plt.savefig(plotfile)
   # show()
    print "Saving to file: %s"%(plotfile)
    stName+=1

    #show()
