import numpy as np
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
C2 = 0.05                     
act = 1
a = (0.01/3600.)*dt;            
b = -1.3 
Pe = 0                        
Ke_larvae = 1
Ke_predator = 1              
attCoeff = 0.18
Nhours =24
Nlarva=1
FishDens = 0.0001	#Predation from fish depends on density of predators
     