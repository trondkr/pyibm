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
NReleaseDatesInit=9500
daysBetweenReleases=7
Nlarva=1
NDaysAlive=30
Nprey=1
initWgt=0.09 #wgt in milligrams
initDepth=25 # Larvae are randomly distributed in range 0-initDepth
randomWgt=1 #1=on, 0=off Initialize weight with random values from initWgt
# Number of nauplii per liter is a function of Nprey time MultiplyPrey: e.g. prey=2*MultiplyPrey
MultiplyPrey=85

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
FishDens = 0.00015	#Predation from fish depends on density of predators
deadThreshold=0.8   #Individuals die if weight is less than 70% of regular weight at length: predation.py (Fiksen et al. 2002)
costRateOfMetabolism=0.5 # The rate of how much full swimming for one time step will cost relative to routine metabolism

"""Here you define how many time steps you want per 24 hours"""
dt_per_day=24
deltaH = 24./(dt_per_day*1.0)	#Hours per timestep

"""Here you define the vertical resolution of behavior and movement meter
deltaZ=0.1 = 10 cm increments, deltaZ=0.05 is 5 cm increments. If resolution is less than 1 meter,
make sure that lastDecimal is larger than zero. lastDecimal=1 : 0.56=>0.6"""
deltaZ=1
lastDecimal=0

"""Generate generic prey size groups based on Daewel et al. 2008 JPR #30 paper.
Here we divide the prey size spectra into 16 groups ranging from 0.1 to 1.6 mm in length.
The width and weight of the prey were derived from the weight and width of Calanus
of equal length. Calculations according to Daewel require length in microgram (converted to
mm at the end for use in ibm)."""
sizeMin=100. # micromm (/1000 to get mm)
sizeMax=1600.
prey_LENGTH=[]
"""Calculate the sizes based on a range of 16 from 100 to 1600 micrograms"""
for i in range(16):
    prey_LENGTH.append(sizeMin*(i+1))

"""Calculate the density function of prey of a given size according to Daewel et al. 2008"""
tm=0 # total biomass
SDpl=np.zeros(16)
SPpl=np.zeros(16)
m=np.zeros(16)

for i in range(16):
    SDpl[i] = 695.73 * np.exp(-0.0083*prey_LENGTH[i])
    m[i] = np.exp( 2.772 * np.log(prey_LENGTH[i]) - 7.476)
    tm = tm+m[i]

"""Here we calulate the percent of each size group represented as a biomass percent"""
for i in range(16):
    SPpl[i] = (SDpl[i] * m[i])/tm

prey_WGT   = [0.1, 0.25, 0.45, 1.0, 1.51, 2.76, 3.76, 6.0, 9.0, 13.24, 15.0, 17.0, 19.0, 23.13, 30.0, 45.0]
prey_WIDTH = [0.08, 0.1, 0.1, 0.1, 0.15, 0.17, 0.18, 0.2, 0.22, 0.25, 0.31, 0.35, 0.38, 0.38, 0.39, 0.40]

prey_LENGTH= np.asarray(prey_LENGTH)/1000. #micromm to mm
prey_D=np.asarray(SPpl)

prey_AREA=np.zeros(len(prey_LENGTH))
"""Calculate the image area (mm^2) of elongated prey"""
for i in range(len(prey_LENGTH)):
    prey_AREA[i]=0.75*(prey_LENGTH[i])*(prey_WIDTH[i])


"""Settings for progressbar"""
empty  =u'\u25FD'
filled =u'\u25FE'
progress  = progressBar(color='red',width=30, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))
