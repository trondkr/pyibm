Module predation

  !
  ! ---------------------------------------------------------------
  ! This module calculates the bioenergetics for larval cod such as metabolism
  !
  ! email : Trond.Kristiansen@imr.no
  ! ---------------------------------------------------------------
  !
  ! USAGE:
  ! mortality, didStarve,StarvationMortality, FishMortality, InvertebrateMortality = predF90.predation.fishpredandstarvation(Mortality,Starved,FishDens,L[cohort,ind,larvaIndex-1]*mm2m,
  !                                                                                              W[cohort,ind,larvaIndex-1],
  !                                                                                              Eb,suming,stomachFullness,beamAttCoeff,m2mm,
  !                                                                                              Ke_predator,deadThreshold,deltaH,dt)
  ! ---------------------------------------------------------------------------------------------------------------------------------------------------

  Implicit None

  ! f2py --verbose --fcompiler=intelem -c -m predF90 perception.f90 predation.f90 --f90flags="-no-heap-arrays"

Contains

  Subroutine FishPredAndStarvation(Mortality,Starved,StarvationMortality, &
  &FishMortality, InvertebrateMortality,FishDens,Larval_m,&
  &Larval_wgt,Eb,suming,stomachFullness,beamAttCoeff,m2mm,&
  &Ke_predator,deadThreshold,dh,seconds)

    USE perception

    double precision Mortality, Starved, FishDens, Larval_m, Larval_wgt, Eb,suming
    double precision beamAttCoeff,m2mm, Ke_predator, deadThreshold,stomachFullness
    double precision Em, VisFieldShape, FishSwimVel, aPred, bPred
    double precision visual, image, EyeSens, pi, LarvalShape, contrast, predator_VISUAL
    double precision StarvationMortality, FishMortality, InvertebrateMortality, LarvalWidth
    Integer setMort, dh, seconds, IER

    !f2py intent(in,out) Mortality, Starved,StarvationMortality, FishMortality, InvertebrateMortality
    !f2py intent(in,overwrite) FishDens,Larval_m,Larval_wgt,Eb,suming
    !f2py intent(in,overwrite) stomachFullness,beamAttCoeff,m2mm,Ke_predator,deadThreshold,dh,seconds


    !This routine takes larval length in meter as input (Larval_m)
    LarvalShape = 0.2		!Larval width:length ratio
    contrast=0.3
    Em = 5.0e4
    VisFieldShape = 0.5
    FishSwimVel = 0.10		!Fish cruising velocity (m/s)
    aPred = 1.78e-5
    bPred = -1.3	!.78d-5=0.1/3600., -1.3 Parameters for purely size-dependent mortality Fiksen et al. 2002
    StarvationMortality=1.0e-6
    setMort=1

    !Values of fish and larval length are all in meters in this subroutine,
    !which differs from the rest of the routines that are all in mm.
    LarvalWidth   = LarvalShape*Larval_m
    image = LarvalWidth*Larval_m*0.75
    pi=4.D0*DATAN(1.D0)

    IER=0
    visual=0.0

    !All input to getr is either in m (or per m), or in mm (or per mm). Here we use meter (m):
    call getr(predator_VISUAL,beamAttCoeff/m2mm,contrast,image,Em,Ke_predator,Eb, IER)

    !Calculate lethal encounter rate with fish setMort is either 0 (off) or 1 (on)
    FishMortality = setMort*(VisFieldShape*pi*(predator_VISUAL**2)*FishSwimVel*FishDens)

    InvertebrateMortality = setMort*OtherPred(Larval_m,aPred,bPred,m2mm)
    Starved = aliveOrDead(Larval_wgt, Larval_m,suming,stomachFullness,deadThreshold,m2mm)
    Mortality = (InvertebrateMortality + FishMortality + Starved*StarvationMortality)*seconds*dh

   ! print *,'inv', InvertebrateMortality*seconds*dh, 'fish',FishMortality*seconds*dh, 'starv',Starved*StarvationMortality*seconds*dh
    return

  end subroutine fishPredAndStarvation

  double precision function OtherPred(L_m,aPred,bPred,m2mm)
    double precision L_m,aPred,bPred,m2mm
    !Calculate death risk from other sources"""
    OtherPred = aPred*(L_m*m2mm)**bPred
    return
  end function OtherPred

  double precision function WeightAtLength(L_m,m2mm)
    double precision L_m,m2mm
    !Calculate the dry body mass (ug) from length in mm (Folkvord 2005)
    WeightAtLength = (exp(-9.38+4.55*log(L_m*m2mm)-0.2046*(log(L_m*m2mm))**2.))*1000.
    return
  end function WeightAtLength

  double precision function aliveOrDead(Larval_wgt, Larval_m,suming,stomachFullness,deadThreshold,m2mm)
    ! For a given length, the weight should be a given reference weight. If the weight is
    ! less than regular weight at length, then we assume staarvation is occurring. If the weight
    ! is less than 70% (deadThreshold - initLarva.py) of the weight it should have at length we assume the larvae is dead
    ! and massively increase the mortality rate to reflect death.
    !
    ! Trond Kristiansen, 02.12.2009, 17.02.2010"""
    double precision Larval_wgt,Larval_m,suming,stomachFullness, refWeight, deadThreshold, m2mm

    refWeight = WeightAtLength(Larval_m,m2mm)

    if (refWeight > Larval_wgt*1000. .OR. ((suming < 0.00000001 .AND. stomachFullness < 0.0000001))) then
        aliveOrDead=1
        !print 'starvation %s refWgt: %s wgt: %s'%(starvation,refWeight, Larval_wgt*1000.)
    else
        aliveOrDead=0.0
    end if

    if (refWeight*deadThreshold > Larval_wgt*1000.) then
        !Give a very high probability of death when belowe 80% (deadThreshold)"""
        aliveOrDead = 10000.0
        !dead = 0
        !print 'Larva died of starvation'
    end if

    return
  end function AliveOrDead

end module predation
