import datetime
import numpy as np
import IOlight
from initLarva import *

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 7, 19)
__modified__ = datetime.datetime(2009, 7, 19)
__version__  = "0.1"
__status__   = "Production"
__notes__    = "Converted from Fortran to Python. Based on growth.f90"

def FishPredAndStarvation(FishDens,Larval_m,Larval_wgt,attCoeff,Ke,Eb,seconds):
    """This routine takes larval length in meter as input (Larval_m)"""
    LarvalShape = 0.2		#Larval width:length ratio
    PreyContrast=0.3
    EyeSens = 5.0e4
    VisFieldShape = 0.5
    FishSwimVel = 0.10		#Fish cruising velocity (m/s)
    aPred = 1.78e-5; bPred = -1.3	#.78d-5=0.1/3600., -1.3 Parameters for purely size-dependent mortality
    Ambush = 0.5					#Fraction of invertebrate predation that is ambush	
    StarvationMortality=1e-6
    setMort=1
       
    """Values of fish and larval length are all in meters in this subroutine,
     which differs from the rest of the routines that are all in mm."""
    
    LarvalWidth = LarvalShape*Larval_m                   # Larval length in meter x larval shape
    #TODO: Check this one from Fiksen et al.
    PreyImageArea = LarvalWidth*Larval_m*0.75		#Larval image area 
    
    """Finding the visual range in m (from AU97)"""
    VisualRange = np.sqrt(EyeSens*PreyContrast*PreyImageArea*(Eb/(Ke+Eb))) #Approximation	     
    
    """Calculate exact visual range (above ca 5 cm)"""
    if VisualRange>0.05:
        VisualRange = IOlight.getPerceptionDistance(EyeSens,attCoeff,Ke,PreyImageArea,Eb)
   
    """Calculate lethal encounter rate with fish setMort is either 0 (off) or 1 (on)"""
    FishMortality = setMort*(VisFieldShape*np.pi*(VisualRange**2)*FishSwimVel*FishDens)
    InvertebrateMortality = (setMort*(Ambush*OtherPred(Larval_m,aPred,bPred) + (1.-Ambush)*OtherPred(Larval_m,aPred,bPred)))
    Starved=aliveOrDead(Larval_wgt, Larval_m)

    Mortality = (InvertebrateMortality + FishMortality + Starved*StarvationMortality)*seconds
    
    #print 'inv', InvertebrateMortality, 'fish',FishMortality, 'starv',aliveOrDead(Larval_wgt, Larval_m)*StarvationMortality
    return Mortality, Starved

def OtherPred(L_m,aPred,bPred):     
    """Calculate death risk from other sources"""
    OtherPred = aPred*(L_m*m2mm)**bPred

    return OtherPred

def WeightAtLength(L_m):
    """Calculate the dry body mass (ug) from length in mm (Folkvord 2005)""" 
    return (np.exp(-9.38+4.55*np.log(L_m*m2mm)-0.2046*(np.log(L_m*m2mm))**2.))*1000.   

def aliveOrDead(Larval_wgt, Larval_m):
    """For a given length, the weight should be a given reference weight. If the weight is
    less than regular weight at length, then we assume staarvation is occurring. If the weight
    is less than 70% (deadThreshold - initLarva.py) of the weight it should have at length we assume the larvae is dead
    and massively increase the mortality rate to reflect death.
    
    Trond Kristiansen, 02.12.2009"""
    
    refWeight = WeightAtLength(Larval_m)
    if (refWeight > Larval_wgt*1000.):
        starvation=1
        #print 'starvation %s refWgt: %s wgt: %s'%(starvation,refWeight, Larval_wgt*1000.)
    else:
        starvation=0
        
    if (refWeight*deadThreshold > Larval_wgt*1000.):
        starvation = 2
        #print 'Larva died of starvation'

    return starvation