import numpy as np
from progressBar import progressBar
"""
=======================
Define global variables
=======================
"""
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
contrast = 0.4 #Inherent contrast after Fiksen,02
mm2m = 0.001
C2 = 0.05                     
act = 1
a = (0.01/3600.)*dt;            
b = -1.3 
Pe = 0                        
Ke_larvae = 1
Ke_predator = 1              
attCoeff = 0.18
Nhours =24
Nlarva=100
FishDens = 0.0001	#Predation from fish depends on density of predators
Nprey=1

""" Initialize Calanus """
calanus_D  = [0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #(#/ltr)
calanus_W  = [0.33, 0.49, 1., 1.51,2.09, 2.76, 4.18, 13.24, 23.13, 63.64, 169.58, 276.29, 276.29] #(micrograms)
calanus_L1 = [0.22, 0.27, 0.4, 0.48, 0.55, 0.61, 0.79, 1.08, 1.38, 1.8, 2.43, 2.11, 2.11] #(length in mm)
calanus_L2 = [0.1, 0.1, 0.1, 0.15,  0.18, 0.2, 0.22, 0.25, 0.31, 0.41, 0.52, 0.65, 0.65] #(width in mm)

"""Settings for progressbar"""
empty  =u'\u25FD'
filled =u'\u25FE'
progress  = progressBar(color='red',width=30, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))

    