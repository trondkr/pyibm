import math

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
        
