#!/usr/bin/env python
from pylab import *
import numpy as np
import string, os, sys
import datetime, types
from netCDF4 import Dataset

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 8, 13)
__modified__ = datetime.datetime(2009, 8, 13)
__version__  = "0.1"
__status__   = "Development"

missingValue=-9.99e-35

station="test_res.nc"
cdf=Dataset(station,"r")

depth=cdf.variables["depth"][:,:,:,:]
time=(cdf.variables["time"][:])
timeIndex=(cdf.variables["timeIndex"][:,:,:,:])

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
ax = fig.add_subplot(3,1,1)
colors=('r','g','b','c','y','m','k','r','g','b','c','y','m','k','r','g','b','c','y','m','k')

"""Calculate 24 hour running mean of specific growth rate"""
Ker = np.ones((24,), dtype='float')
Ker[0], Ker[-1] = 0.5, 0.5
Ker = Ker/(24.)
depth24H=np.zeros(depth.shape,dtype=float64)

Ncohorts=11#len(depth[:,0,0,0])-4
Nindividuals=len(depth[0,:,0,0])
print "Number of individuals %s and cohorts %s"%(Nindividuals, Ncohorts)
for cohort in range(Ncohorts):
    index1=int(timeIndex[cohort,0,0,0])+1
    index2=int(timeIndex[cohort,0,1,0])-1
    
    days=int((index2-index1)/24)
    print "Finding %s days period of growths for cohort %s"%(days,cohort)
    
    meanSGR24H=np.zeros((days),dtype=float64)
    meanIndSGR24H=np.zeros((Nindividuals,days),dtype=float64)
    
    stdSGR24H=np.zeros((days),dtype=float64)
     
    time24H=np.zeros((days),dtype=float64)
    
    for i in range(days):
        k1=index1+i*24
        k2=index1+(i+1)*24
       
        meanSGR24H[i]=sum(sgr[cohort,:,k1:k2,0])/Nindividuals
     
        j=0
        
        for j in range(Nindividuals):
            meanIndSGR24H[j,i]=sum(sgr[cohort,j,k1:k2,0])
        
        stdSGR24H[i]=std(meanIndSGR24H[:,i])
        time24H[i]=time[k1]
    
    stdMinus = meanSGR24H - stdSGR24H
    stdPlus  = meanSGR24H + stdSGR24H
    
 #   sgr24H[cohort,ind,index1:index2,0]= np.convolve(sgr[cohort,ind,index1:index2,0], Ker, mode='same')
    ax.plot(time24H,meanSGR24H,colors[cohort],linewidth = 3)
    ax.fill_between(time24H, stdMinus, stdPlus, facecolor='black',alpha=0.5)

    #ax.plot(time[index1:index2],sgr[cohort,ind,index1:index2,0],colors[cohort],linewidth = 2)
    
ax.set_title('Specific growth rate')


ax = fig.add_subplot(3,1,2)
for cohort in range(Ncohorts):

     index1=int(timeIndex[cohort,0,0,0])+1
     index2=int(timeIndex[cohort,0,1,0])-1
     meanWGT = mean(wgt[cohort,:,index1:index2,0],0)
     
     ax.plot(time[index1:index2],meanWGT,colors[cohort],linewidth = 2)
     stdMinus=mean(wgt[cohort,:,index1:index2,0],0) - std(wgt[cohort,:,index1:index2,0],0)
     stdPlus=mean(wgt[cohort,:,index1:index2,0],0) + std(wgt[cohort,:,index1:index2,0],0)
     ax.fill_between(time[index1:index2], stdMinus, stdPlus, facecolor='black',alpha=0.5)
     
ax.set_title('Weight')


ax = fig.add_subplot(3,1,3)
for cohort in range(Ncohorts):#len(depth[:,0,0,0])-1)
        
    index1=int(timeIndex[cohort,0,0,0])+1
    index2=int(timeIndex[cohort,0,1,0])-1
    meanDEPTH = mean(-depth[cohort,:,index1:index2,0],0)
 
    ax.plot(time[index1:index2],meanDEPTH,colors[cohort],linewidth = 2)
    stdMinus=mean(-depth[cohort,:,index1:index2,0],0) - std(-depth[cohort,:,index1:index2,0],0)
    stdPlus=mean(-depth[cohort,:,index1:index2,0],0) + std(-depth[cohort,:,index1:index2,0],0)
    ax.fill_between(time[index1:index2], stdMinus, stdPlus, facecolor='black',alpha=0.5)
    
ax.set_title('Weight')


show()