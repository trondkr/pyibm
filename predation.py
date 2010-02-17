import datetime
import numpy as np
#import IOlight
from initLarva import *
import calclight
import perception

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2009, 7, 19)
__modified__ = datetime.datetime(2010, 2, 17)
__version__  = "0.1"
__status__   = "Production"
__notes__    = "Converted from Fortran to Python. Based on growth.f90. Edited 19.07.2009, 17.02.2010"

def FishPredAndStarvation(grdSTATION,dh,FishDens,Larval_m,Larval_wgt,attCoeff,Eb,seconds,ing,stomachFullness):
    """This routine takes larval length in meter as input (Larval_m)"""
    LarvalShape = 0.2		#Larval width:length ratio
    PreyContrast=0.3
    EyeSens = 5.0e4
    VisFieldShape = 0.5
    FishSwimVel = 0.10		#Fish cruising velocity (m/s)
    aPred = 1.78e-5; bPred = -1.3	#.78d-5=0.1/3600., -1.3 Parameters for purely size-dependent mortality Fiksen et al. 2002
    StarvationMortality=1.0e-6
    setMort=1
       
    """Values of fish and larval length are all in meters in this subroutine,
     which differs from the rest of the routines that are all in mm."""
    LarvalWidth   = LarvalShape*Larval_m
    PreyImageArea = LarvalWidth*Larval_m*0.75	
   
    IER=0; visualRange=0.0 
    """All input to getr is either in m (or per m), or in mm (or per mm)"""
    visualRange, IER = perception.perception.getr(visualRange,beamAttCoeff/m2mm,PreyContrast,PreyImageArea,EyeSens,Ke_predator,Eb, IER)
    
    """Calculate lethal encounter rate with fish setMort is either 0 (off) or 1 (on)"""
    FishMortality = setMort*(VisFieldShape*np.pi*(visualRange**2)*FishSwimVel*FishDens)
    InvertebrateMortality = setMort*OtherPred(Larval_m,aPred,bPred)
    Starved, dead = aliveOrDead(Larval_wgt, Larval_m,ing,stomachFullness)

    Mortality = (InvertebrateMortality + FishMortality + Starved*StarvationMortality)*seconds*dh
    
    #print 'inv', InvertebrateMortality, 'fish',FishMortality, 'starv',Starved*StarvationMortality
    return Mortality, Starved, dead

def OtherPred(L_m,aPred,bPred):     
    """Calculate death risk from other sources"""
    OtherPred = aPred*(L_m*m2mm)**bPred
    return OtherPred
    
def WeightAtLength(L_m):
    """Calculate the dry body mass (ug) from length in mm (Folkvord 2005)""" 
    return (np.exp(-9.38+4.55*np.log(L_m*m2mm)-0.2046*(np.log(L_m*m2mm))**2.))*1000.   

def aliveOrDead(Larval_wgt, Larval_m,ing,stomachFullness):
    """For a given length, the weight should be a given reference weight. If the weight is
    less than regular weight at length, then we assume staarvation is occurring. If the weight
    is less than 70% (deadThreshold - initLarva.py) of the weight it should have at length we assume the larvae is dead
    and massively increase the mortality rate to reflect death.
    
    Trond Kristiansen, 02.12.2009, 17.02.2010"""
    
    refWeight = WeightAtLength(Larval_m)
    dead=1
    if (refWeight > Larval_wgt*1000. or (sum(ing) < 0.00000001 and stomachFullness < 0.0000001)):
        starvation=1
        #print 'starvation %s refWgt: %s wgt: %s'%(starvation,refWeight, Larval_wgt*1000.)
    else:
        starvation=0
        
    if (refWeight*deadThreshold > Larval_wgt*1000.):
        """Give a very high probability of death when belowe 80% (deadThreshold)"""
        starvation = 100000
        dead = 0
        #print 'Larva died of starvation'

    return starvation, dead