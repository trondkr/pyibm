import math, datetime

import IOtime

"""
Functions to estimate maximum surface light as a function of day of the year and latitude. 
Based on the surface light and depth, and size of the larvae we calculate the perception distance of fish.

Converted from Fortran to Python on 10.06.2008 by Trond Kristiansen
Trond.Kristiansen@imr.no

Functions in module:

1. [umol s-1 m-2]=isurface_light(julian, Lat, hour)
2. [mm]=get_perception_distance(k,Ke,Ap,Eb)

"""  
__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 6, 10)
__modified__ = datetime.datetime(2008, 6, 10)
__version__  = "1.1"
__status__   = "Production"

def surface_light(julian, Lat, hour):
    
    """
    Max light at sea surface
    """
    MAXLIG = 500
    
    # Need day and hour of year
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    #gtime = gregorian(julian)
    gtime = IOtime.julian_to_date(julian,hour)
    
    D = float(gtime[2] + sum(days_in_month[0:int(gtime[1])-1]))
    H = float(gtime[3])
    
    P = math.pi
    TWLIGHT = 5.76
    
    MAXLIG = MAXLIG * 0.1 + MAXLIG * abs(math.sin(math.pi*(D*1.0)/365.0))
    
    """
    Originally the calculation of light as a function of lat, time of day, and time of year was setup for
    sind, cosd of degrees. I converted all to radians using math.radians to suit the python standard.
       
    Test suite: 2440000 should give 1968, 05, 23, 00, 00, 00
    for k in range(24):
        surface_light(2440000+k/24.,42.0)
        print 2440000+k/24.
        
    print julian_day(1968,05,23,00)
    """

    DELTA = 0.3979*math.sin(math.radians(0.9856*(D-80)+ 1.9171*(math.sin(math.radians(0.9856*D))-0.98112)))
    H12 = DELTA*math.sin(math.radians(Lat*1.))- math.sqrt(1.-DELTA**2)*math.cos(math.radians(Lat*1.))*math.cos(math.radians(15.0*12.0))
    HEIGHT = DELTA*math.sin(math.radians(Lat*1.))- math.sqrt(1.-DELTA**2)*math.cos(math.radians(Lat*1.))*math.cos(math.radians(15.0*H))
      
    V = math.asin(HEIGHT)
    
    if (V >= 0):                 
        s_light = MAXLIG*(HEIGHT/H12) + TWLIGHT
    elif (V >= -6):
        s_light = ((TWLIGHT - 0.048)/6.)*(6.+V)+.048
    elif (V >= -12):
        s_light = ((0.048 - 1.15e-4)/6.)*(12.+V)+1.15e-4
    elif (V >= -18):
        s_light = (((1.15e-4)-1.15e-5)/6.)*(18.+V)+1.15e-5
    else:
        s_light = 1.15e-5
        
    return s_light
    
def get_perception_distance(k,Ke,Ap,Eb):
    """
    Converted getr.f90 from Fortran to python. Based on Dag Aksnes formula for
    iterating over perception distance.
    Trond Kristiansen 10.06.2008
    Trond.Kristiansen@imr.no
    """
    c = 3*k
    C0 = 0.4 #Inherent contrast after Fiksen,02
    mm2m = 0.001
    Vc = 10000 #Size-specific sensitivity of the visual system (Fiksen,02)
    
    #Initial guess of visual range (RST)
    R2= abs(C0)*Ap*Vc*(Eb/(Ke+Eb))
    RST = math.sqrt(R2)
    
    #Upper boundary of allowed error of visual range
    EPS = .0000001
    #Maximum number of iteration steps
    IEND = 100
    
    #Prepare iteration
    r = RST
    TOL = r
    FR2=math.log(abs(C0)*Ap*Vc)
    FR1=math.log(((Ke+Eb)/Eb)*r*r*math.exp(c*r))
    F1 = FR1-FR2
    FDER = c + 2./r 
    TOLF = 100. * EPS
    
    for I in range(IEND):
        if (abs(F1)>0): # if 1
            if (abs(FDER)>0): # if 2
                DX = F1/FDER
                r = r - DX
                if (r>=0): # if 3
                    TOL = r
                    FR2=math.log(abs(C0)*Ap*Vc)
                    FR1=math.log(((Ke+Eb)/Eb)*r*r*math.exp(c*r))
                    F1 = FR1-FR2
                    FDER = c + 2./r 
                    TOL = EPS
                    AS = abs(r)
                    if (AS-1>0): # if: 4
                        TOL = TOL*AS;
                  
                    if (abs(DX)-TOL<=0): # if: 5
                        if (abs(F1)-TOLF>0):   # if: 6
                            continue
                        else: # if: 6
                            break
                    
                    else: # if: 5
                        break
                  
                else: # if: 3
                    break
                    
            else: # if: 2
                break
        else: # if: 1
            break
        
    return r/mm2m # returns mm
        

    
   