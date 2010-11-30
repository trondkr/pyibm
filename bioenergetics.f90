Module bioenergetics

  !
  ! ---------------------------------------------------------------
  ! This module calculates the bioenergetics for larval cod such as metabolism
  ! ingestion, feeding, growth rates.
  !
  ! email : Trond.Kristiansen@imr.no
  ! ---------------------------------------------------------------

  Implicit None
  
  ! f2py --verbose --fcompiler=intelem -c -m bioenergetics perception.f90 bioenergetics.f90        
Contains
  

  
  Subroutine growth (Tdata,NLGdata,windX,windY,&
  &Spre,prey_AREA,prey_LENGTH,prey_D,&
  &prey_WGT,Larval_mm,Larval_wgt,sec2day,&
  &mg2ug,ltr2mm3, m2mm,mm2m,&
  &travel,tau,micro2m,Eb,&
  &contrast,Ke_larvae,beamAttCoeff,deltaH,&
  &dt,II,depth,prey,gut_size)

    USE perception

    Real meta, Larval_mm, Tdata, Larval_wgt, GR, GR_gram, assi, stomachFullness,Spre
    Real sec2day, mg2ug, ltr2mm3, m2mm, mm2m, travel, tau, micro2m, gut_size, pi, larval_VISUAL, image
    Real zoop, NLGdata
    Real Eb, Em, contrast, Ke_larvae, beamAttCoeff, windX, windY, omega, eps
    Integer deltaH, dt, II, depth, prey, j, IER
    
    double precision, dimension(II) :: ing, enc, hand, pca, visual, prey_AREA, prey_LENGTH, prey_D, prey_WGT
                         
!f2py intent(in) Tdata,NLGdata,windX,windY,Spre,prey_AREA,prey_LENGTH,prey_D,prey_WGT,Larval_mm,Larval_wgt,sec2day,mg2ug,ltr2mm3,m2mm,mm2m,travel,tau,micro2m,Eb,contrast,Ke_larvae,beamAttCoeff,deltaH,dt,II,depth,prey,gut_size
!f2py intent(out) ing,GR_gram,GR,meta,assi,stomachFullness,zoop
    print *,Tdata,NLGdata,windX,windY
    print *,Spre,prey_AREA,prey_LENGTH
    print *,prey_D,prey_WGT,Larval_mm,Larval_wgt
    print *,sec2day, mg2ug,ltr2mm3, m2mm,mm2m
    print *,travel,tau,micro2m,Eb,contrast,Ke_larvae
    print *,beamAttCoeff,deltaH,dt,II,depth,prey,gut_size
    
    ! Mouthsize of larvae. The larvae can only capture calanus of sizes less
    ! than the larval mouthsize. Folkvord et al.
    !mouth = exp(-3.27+1.818*log(Larval_mm)-0.1219*(log(Larval_mm))**2.)
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
    ! Calculate growth rate
    GR = max(0.0, sec2day*log(0.01*(1.08 + 1.79*Tdata - 0.074*Tdata*log(Larval_wgt) &
    &- 0.0965*Tdata*log(Larval_wgt)**2 &
    &+ 0.0112*Tdata*log(Larval_wgt)**3) + 1)*dt*deltaH)
    
    ! Growth rate converted to percentage
    GR_gram = (exp(GR) - 1.0)*Larval_wgt


    ! Use large phytoplankton as proxy for zooplankton
    zoop=0.45 * NLGdata * 1.e4

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
        visual(j)=sqrt(Em*contrast*(prey_AREA(j))*(Eb/(Ke_larvae+Eb)))
        ! All input to getr is either in m (or per m), or in mm (or per mm)"""
        image=prey_AREA(j)
        call getr(larval_VISUAL,beamAttCoeff/m2mm,contrast,image,Em,Ke_larvae,Eb, IER)
        visual(j)=larval_VISUAL
        pca(j) = max(0.0,min(1.0,-16.7*(prey_LENGTH(j)/Larval_mm) + 3.0/2.0))
        
        omega = 1.9*(eps*visual(j)*mm2m)**0.667
        omega =  omega * m2mm ! From m/s to mm/s
    
        ! Calculate handling time, encounter rate, and probability of capture"""
        hand(j) = 0.264*10**(7.0151*(prey_LENGTH(j)/Larval_mm)) ! Walton 1992
        enc(j) = ((0.667*pi*(visual(j)**3.)*travel + pi*(visual(j)**2.)*sqrt(prey_LENGTH(j)**2.+ 2.*omega**2.)*travel*tau) &
        &* (prey_D(j)*((prey+1)*zoop))* ltr2mm3)
        ing(j) = (dt*enc(j)*pca(j)*prey_WGT(j)*micro2m / (1 + hand(j)))*deltaH

        print *,"j:",j, " enc:",enc(j)*dt*deltaH, " mm:",Larval_mm, "pca:",pca(j), "SPpl:",prey_D(j)
        print *,"j:",j, " ing:",ing(j), "visual:",visual(j)
    end do
        
    print *,"prey density: ", sum(prey_D(:)*((prey+1)*zoop))
  
    ! Calculate stomach fullness
    stomachFullness =  (min(gut_size*Larval_wgt,Spre + sum(ing)))/(Larval_wgt*gut_size)
    return 
    
    end subroutine growth
end module bioenergetics