import numpy as np
from progressBar import progressBar
import datetime as datetime

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 9)
__modified__ = datetime.datetime(2009, 8, 10)
__version__  = "1.1.5"
__status__   = "Production"

"""
====================================================================
Define global 
Variables that are imported to ibm.py using from initLarva import *
====================================================================
"""

"""Number of release dates and cohorts. This will be lowered if NReleaseDatesInit*daysBetweenReleases is
more than total number of simulation days : see function init in ibm.py"""
NReleaseDatesInit=100
daysBetweenReleases=30
Nlarva=1
NDaysAlive=30
Nprey=3
initWgt=0.090 #wgt in milligrams
initDepth=10
randomWgt=1 #1=on, 0=off Initialize weight with random values from initWgt
# Number of nauplii per liter is a function of Nprey time MultiplyPrey: e.g. prey=2*MultiplyPrey
MultiplyPrey=10

missingValue=-9.99e-35
dt      = 3600                   
sec2day = 1.0/86400.0
gut_size= 0.06
stomach_threshold = 0.3
tau   = 2.0                     
omega = 0.0
f     = 0.43                 
pi    = np.pi
mm2m  = 0.001
m2mm  = 1000.
ltr2mm3 = 1e-6
micro2m = 0.001
contrast = 0.3 #Inherent contrast after Fiksen,02
mm2m = 0.001
mg2ug=1000.0
C2 = 0.05                     
act = 1
a = (0.01/3600.)*dt;            
b = -1.3 
Pe = 0                        
Ke_larvae = 1
Ke_predator = 1              
attCoeff = 0.18
beamAttCoeff=attCoeff*3.0
Nhours =24
FishDens = 0.0001	#Predation from fish depends on density of predators
deadThreshold=0.8   #Individuals die if weight is less than 70% of regular weight at length: predation.py

"""Here you define how many time steps you want per 24 hours"""
dt_per_day=48
deltaH = 24./(dt_per_day*1.0)	#Hours per timestep
   
"""Here you define the vertical resolution of behavior and movement meter
deltaZ=0.1 = 10 cm increments, deltaZ=0.05 is 5 cm increments. If resolution is less than 1 meter,
make sure that lastDecimal is larger than zero. lastDecimal=1 : 0.56=>0.6"""
deltaZ=1.0
lastDecimal=0

""" Initialize Calanus """
calanus_D  = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #(#/ltr)
calanus_W  = [0.33, 0.49, 1., 1.51,2.09, 2.76, 4.18, 13.24, 23.13, 63.64, 169.58, 276.29, 276.29] #(micrograms)
calanus_L1 = [0.22, 0.27, 0.4, 0.48, 0.55, 0.61, 0.79, 1.08, 1.38, 1.8, 2.43, 2.11, 2.11] #(length in mm)
calanus_L2 = [0.1, 0.1, 0.1, 0.15,  0.18, 0.2, 0.22, 0.25, 0.31, 0.41, 0.52, 0.65, 0.65] #(width in mm)
calanus_Area=np.zeros(len(calanus_L2))
"""Calculate the image area (mm^2) of elongated prey"""
for i in range(len(calanus_L1)):
    calanus_Area[i]=0.75*(calanus_L1[i])*(calanus_L2[i])
    
"""Settings for progressbar"""
empty  =u'\u25FD'
filled =u'\u25FE'
progress  = progressBar(color='red',width=30, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))

    
