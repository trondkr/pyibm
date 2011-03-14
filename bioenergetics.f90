Module bioenergetics

  !
  ! ---------------------------------------------------------------
  ! This module calculates the bioenergetics for larval cod such as metabolism
  ! ingestion, feeding, growth rates. The return values are:
  ! suming,ing,GR_gram,meta,assi,stomachFullness,zoop where
  ! suming = sum(ing(:))  - sum function faster to do in Fortran than Python
  ! ing(0:II) - an array contaiing the ingestin for each prey size 0 to II
  ! GR_gram - growth rate in weight
  ! GR growth rate
  ! meta - metabolism
  ! assi - assimilation
  ! stomachFullness - stomach fullness
  ! zoop - zooplankton density
  !
  ! email : tk@trondkristiansen.com
  ! ---------------------------------------------------------------
  !
  ! USAGE:
  !   suming,ing,GR_mg,meta,assi,stomachFullness,zoop = bioenergetics.bioenergetics.growth(suming,ing,GR_mg,meta,assi,stomachFullness,zoop,
  !                                                                                      Tdata,windX,windY,
  !                                                                                      S[cohort,ind,larvaIndex-1],prey_AREA,prey_LENGTH,
  !                                                                                      prey_D,prey_WGT,L[cohort,ind,larvaIndex-1],W[cohort,ind,larvaIndex-1],
  !                                                                                      sec2day, mg2ug,ltr2mm3, m2mm,mm2m,f,tau,ug2mg,Eb,contrast,Ke_larvae,
  !                                                                                      beamAttCoeff,MultiplyPrey,deltaH,dt,depth,gut_size,II)
  ! ---------------------------------------------------------------------------------------------------------------------------------------------------

  Implicit None

  ! f2py --verbose --fcompiler=intelem -c -m bioenergetics perception.f90 bioenergetics.f90 --f90flags="-no-heap-arrays"
  ! Also in .bash_profile it pays off substantially to only use :export FFLAGS="-m64 -O3 -fast -arch x86_64

Contains

  Subroutine growth (suming,ing,GR_mg,meta,assi,stomachFullness,zoop,&
  &Tdata,windX,windY,&
  &Spre,prey_AREA,prey_LENGTH,prey_D,&
  &prey_WGT,Larval_mm,Larval_wgt,sec2day,&
  &mg2ug,ltr2mm3, m2mm,mm2m,&
  &travel,tau,ug2mg,Eb,&
  &contrast,Ke_larvae,beamAttCoeff,MultiplyPrey,deltaH,&
  &dt,depth,gut_size,II)

    USE perception

    Real meta, Larval_mm, Tdata, Larval_wgt, suming, GR, GR_mg, assi, stomachFullness,Spre
    Real sec2day, mg2ug, ltr2mm3, m2mm, mm2m, travel, tau, ug2mg, gut_size, pi, larval_VISUAL, image
    Real MultiplyPrey
    Real Eb, Em, contrast, Ke_larvae, beamAttCoeff, windX, windY, omega, eps, g
    Integer deltaH, dt, depth, II, j, IER

    double precision, dimension(II) :: ing, enc, hand, pca, visual, prey_AREA, prey_LENGTH, prey_D, prey_WGT, zoop

!f2py intent(in) ing,Tdata,windX,windY
!f2py intent(in) Spre,prey_AREA,prey_LENGTH,prey_D,prey_WGT,Larval_mm,Larval_wgt,sec2day,mg2ug
!f2py intent(in) ltr2mm3,m2mm,mm2m,travel,tau,ug2mg,Eb,contrast,Ke_larvae
!f2py intent(in) beamAttCoeff,MultiplyPrey,deltaH,dt,depth,gut_size,II,zoop
!f2py intent(in,out,overwrite) suming,ing,GR_mg,meta,assi,stomachFullness

    ! NOTE: input units to this subroutine is larval_weight in mg. Output
    ! converted to microgram.

    ! Calculate metabolism
    meta = dt*2.38e-7*exp(0.088*Tdata)*((Larval_wgt*mg2ug)**(0.9)*0.001)*deltaH
    ! Increase metabolism during active hours (light above threshold) Lough et  al. 2005
    if (Eb > 0.001) then
        if (Larval_mm > 5.5) then
            meta = (2.5*meta)
        else
            meta = (1.4*meta)
        end if
    end if
    !Calculate assimilation efficiency
    assi=0.8*(1-0.400*exp(-0.002*(Larval_wgt*mg2ug-50.0)))
    ! Calculate daily growth rate (SGR in percent %)
    GR = 1.08 + 1.79*Tdata - 0.074*Tdata*log(Larval_wgt) &
        &- 0.0965*Tdata*log(Larval_wgt)**2 &
        &+ 0.0112*Tdata*log(Larval_wgt)**3

    ! Growth rate (g) converted to milligram weight (GR_mg) per timestep:
    g =  (alog(GR/100.+1))*sec2day*dt*deltaH
    GR_mg= (Larval_wgt*(Exp(g)-1.))

    ! VISUAL == PERCEPTION of PREY calculations
    Em = (Larval_mm**2.0)/(contrast*0.1*0.2*0.75) !Size-specific sensitivity of the visual system (Fiksen,02)
    ing(:)  =0.0
    enc(:)  =0.0
    hand(:) =0.0
    pca(:)  =0.0
    pi=4.D0*DATAN(1.D0)

    ! Calculate turbulence based on wind stress"""
    eps=(5.82*1.E-9*((sqrt(windX**2 + windY**2)))**3.)/(depth+0.1)

    do j=1,II
        IER=0
        larval_VISUAL=sqrt(Em*contrast*(prey_AREA(j))*(Eb/(Ke_larvae+Eb)))
        ! All input to getr is either in m (or per m), or in mm (or per mm)"""
        image=prey_AREA(j)
        call getr(larval_VISUAL,beamAttCoeff/m2mm,contrast,image,Em,Ke_larvae,Eb, IER)
        visual(j)=larval_VISUAL
        pca(j) = max(0.0,min(1.0,-16.7*(prey_LENGTH(j)/Larval_mm) + 3.0/2.0))

        omega = 1.9*(eps*visual(j)*mm2m)**0.667
        omega =  omega * m2mm ! From m/s to mm/s

        ! Calculate handling time, encounter rate, and probability of capture"""
        hand(j) = 0.264*10**(7.0151*(prey_LENGTH(j)/Larval_mm)) ! Walton 1992
        enc(j) = ((0.667*pi*(visual(j)**3.)*travel + pi*(visual(j)**2.)*sqrt(prey_LENGTH(j)**2.+ 2.*omega**2.)*travel*tau)*(MultiplyPrey*zoop(j))* ltr2mm3)
        ing(j) = (dt*enc(j)*pca(j)*prey_WGT(j)*ug2mg / (1 + hand(j)))*deltaH

        !print *,"j:",j, " enc:",enc(j)*dt*deltaH, " mm:",Larval_mm, "pca:",pca(j), "SPpl:",prey_D(j)
        !print *,"j:",j, " ing:",ing(j), "suming:",sum(ing), "visual:",larval_VISUAL
    end do

    ! Calculate stomach fullness
    stomachFullness =  (min(gut_size*Larval_wgt,Spre + sum(ing)))/(Larval_wgt*gut_size)
    suming=sum(ing)

    return

    end subroutine growth

end module bioenergetics