import datetime
import numpy as np
import IOlight

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 7, 19)
__modified__ = datetime.datetime(2009, 7, 19)
__version__  = "0.1"
__status__   = "Production"
__notes__    = "Converted from Fortran to Python. Based on growth.f90"

def FishPred(FishDens,Larval_m,attCoeff,Ke,Eb,seconds):
    """This routine takes larval length in meter as input (Larval_m)"""
    LarvalShape = 0.2		#Larval width:length ratio
    PreyContrast=0.3
    EyeSens = 5.0e4
    VisFieldShape = 0.5
    FishSwimVel = 0.10		#Fish cruising velocity (m/s)
    aPred = 1.78e-5; bPred = -1.3	#.78d-5=0.1/3600., -1.3 Parameters for purely size-dependent mortality
    Ambush = 0.5					#Fraction of invertebrate predation that is ambush	
    Starvation_mortality=1e-6
    setMort=1
       
    """Values of fish and larval length are all in meters in this subroutine,
     which differs from the rest of the routines that are all in mm."""
    
    LarvalWidth = LarvalShape*Larval_m                   # Larval length in meter x larval shape
    #TODO: Check this one from Fiksen et al.
    PreyImageArea = LarvalWidth*Larval_m*0.75		#Larval image area 
    
    """Finding the visual range in m (from AU97)"""
    VisualRange = np.sqrt(EyeSens*PreyContrast*PreyImageArea*(Eb/(Ke+Eb))) #Approximation	     
    print 'visual range in cm',VisualRange*100.      
    """ Calculate exact visual range (above ca 5 cm)"""
    if VisualRange>0.05:	
       VisualRange = IOlight.getPerceptionDistance(attCoeff,Ke,PreyImageArea,Eb)
    print 'visual range in cm',VisualRange*100. 
    """Calculate lethal encounter rate with fish
     setMort is either 0 (off) or 1 (on)"""
    FishMortality = setMort*(VisFieldShape*np.pi*(VisualRange**2)*FishSwimVel*FishDens)*seconds
    InvertebrateMortality = setMort*(Ambush*OtherPred(Larval_m,aPred,bPred) + (1.-Ambush)*OtherPred(Larval_m,aPred,bPred))

    #TODO fix starvation mortality
    Mortality = (InvertebrateMortality + FishMortality ) #+ ((this%starvation)*Starvation_mortality))*seconds
    #print Mortality
    return Mortality

    
"""Calculate death risk from other sources"""
def OtherPred(L_m,aPred,bPred): 

    OtherPred = aPred*(L_m*1e3)**bPred

    return OtherPred