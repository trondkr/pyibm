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
  !   suming,ing,GR_mg,meta,assi,stomachFullness,zoop = bioenergetics.bioenergetics.growth(pca,psa,suming,ing,GR_mg,meta,assi,stomachFullness,zoop,
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

    double precision meta, Larval_mm, Tdata, Larval_wgt, suming, GR, GR_mg, assi, stomachFullness,Spre
    double precision sec2day, mg2ug, ltr2mm3, m2mm, mm2m, travel, tau, ug2mg, gut_size, pi, visual, image
    double precision MultiplyPrey
    double precision Eb, Em, contrast, Ke_larvae, beamAttCoeff, windX, windY, omega, eps, g
    Integer deltaH, dt, depth, II, j, IER
    double precision m_teta, var_teta, x_star, teta, capture, r_rand
    double precision speed_fish, speed_prey, w, va,f, d, i, t, ats, N
    double precision gape, dt_pca, dt_num, c, m, d_crit, rs, v, capt_psa, capt_pca, dr, pt

    double precision, dimension(II) :: ing, enc, hand, prey_AREA, prey_LENGTH, prey_D, prey_WGT, zoop, pca, psa

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
    !g =  (log((GR+0.0000000001)/100.+1))*sec2day*dt*deltaH
    !print*,"correct: does this fix it",g, meta, assi, Tdata, Larval_wgt
    ! TODO:
    ! In some rare cases with very cold water which is outside the limnits
    ! of this function, I got values for GR less than 1. Taking the log of less than 1
    ! or 0 gives inf (nan) values. I added a check (max(0.0,GR)) to avoid this
    ! problem. But this should be fixed...
    g =  max(0.0,(dlog(GR/100.+1))*sec2day*dt*deltaH)

    !print*,"wrong: should be nan",g,GR,(log((GR+0.0000000001)/100.+1))*sec2day*dt*deltaH,meta, assi, Tdata, Larval_wgt

    GR_mg= (Larval_wgt*(Exp(g)-1.))

    ! VISUAL == PERCEPTION of PREY calculations
    Em = (Larval_mm**2.0)/(contrast*0.1*0.2*0.75) !Size-specific sensitivity of the visual system (Fiksen,02)
    ing(:)  =0.0
    enc(:)  =0.0
    hand(:) =0.0
    pca(:)  =0.0
    pi=4.D0*DATAN(1.D0)

    dt_num=100
    dt_pca=0.1
    dr=0.0
    pt=0.0
    m_teta = pi / 6. !average escape direction
    var_teta = pi / 6. !m_teta*var_teta = standard deviation of escape direction
    x_star = 0.5 * Larval_mm
    ats=0
    speed_fish = 10.0 ! larva/juvenile swimming speed*BodyLength
    speed_prey = 100.0 ! prey swimming speed*BodyLength
    va = speed_fish * Larval_mm

    ! Probabilites are estimated from 100 capture and approach trials
    ! Fiksen and MacKenzie 2002
    N=10

    ! Calculate turbulence based on wind stress"""
    eps=(5.82*1.E-9*((sqrt(windX**2 + windY**2)))**3.)/(depth+0.1)

    ! The function for gape is for cod, but was compared to Brodeur 1997
    ! which looked at diet of walleye pollock. The results of Brodeur are almost perfect with
    ! the estimates from this function atleast for 65 mm pollock.
    ! Brodeur 1997 Environmental Biology of Fishes 51
    gape=exp(-3.720 + 1.818*log(Larval_mm)-0.1219*(log(Larval_mm))**2) !0.128*Larval_mm**(0.923)

    do j=1,II
        IER=0
        visual=sqrt(Em*contrast*(prey_AREA(j))*(Eb/(Ke_larvae+Eb)))
        ! All input to getr is either in m (or per m), or in mm (or per mm)"""
        image=prey_AREA(j)
        call getr(visual,beamAttCoeff/m2mm,contrast,image,Em,Ke_larvae,Eb, IER)

        pca(j) = max(0.0,min(1.0,-16.7*(prey_LENGTH(j)/Larval_mm) + 3.0/2.0))


        ! PCA and PSA ---------------------------------------------------------
        ! Capture and approach probabilites following Fiksen and MacKenzie 2002
        c = 0.5 * gape !Head radie = 0.5* mouthsize
        rs = c+0.1*Larval_mm !if rs = 0.2*r=0.2*0.8*mm (MK2000), NB! rs-c is strike distance!
        d_crit = 0.264 / prey_LENGTH(j) !Ki¿rboe 1999 (Acartia) d_crit :deformation rate (s-1) at which prey detects the predator
        w = speed_prey * prey_LENGTH(j)
        capt_pca=0.0
        capt_psa=0.0
        travel=0.43
        pt=0.0
        ats=0

        ! Calculate the probability of approach and capture
        enc_loop: Do i = 1, N !Iterates over N encounters
            d = max(visual, c)!Initial predator-prey distance
            approach_loop: Do t = 1, dt_num !Approach from r to rs

               if (d > rs) Then !Maximum non-detective velocity
                  v = ((d_crit*2.*(d**4.))/(3.*c*((d**2.)-(c**2.))))
                  v = dt_pca * min (1.*Larval_mm, v)
               else
                  capt_psa = capt_psa + 1.
                  exit approach_loop
               end if
               !
               dr = - v * (1.-(3*c/(2.*d))+(c**2)/(2*d**3))!Rate of approach (mm/dt)
               d = d + dr !Approach in dt (mm)
               !
            end do approach_loop
            !
            pt = pt + t*dt_pca !Minimum pursuit time
            !
            ! Create a relationship between lp/mm, tp and hour of day
            !
     110    Call random_number (r_rand)
            !
            var_teta = n_dev (r_rand)

            ! Escape direction of prey is 30 deg +- 30 deg (pi/6+-pi/6)
            ! Use the normal distribution (-1 0 1), multiply with pi/6
            ! and add pi/6 to get the correct distribution with mean
            ! value of 0.52 (pi/6).
            !
     666    teta = (m_teta -  var_teta * m_teta)

            If (teta > pi) teta = 2.0 * pi - teta
              teta = Abs (teta)!Axially homogenous
            !
            ats = ats + 1 !Attempts
            !
            if (ats > 3) Go To 333
            if (teta < pi*0.5) Then
               if (gape*0.5/x_star > tan(teta)) Then
                  if ((x_star*cos(teta)/w) < ((rs-c+x_star)/va)) Then
                     Go To 110
                  end if
               end if
            end if
            !
            capture = (w/va) * (sin(teta)*(rs+c)+(gape/2.)*cos(teta))!Capture model
            if (capture < gape*0.5) Then
               capt_pca = capt_pca + 1. !Prey inside the area swept by the larvae
               !print*,"CAPT PCA",capt_pca
            else
               Go To 110
            end If
     333    ats = 0
            !
         end Do enc_loop
        ! Results of iterations
        pt = pt / N * 1.

        psa(j) = min(1.0, capt_psa / N * 1.) !N*1.  !Successful approaches
        pca(j) = min(1.0, capt_pca / N * 1.) !N*1.	!Successful attacks

        ! PCA done -------
        omega = 1.9*(eps*visual*mm2m)**0.667
        omega =  omega * m2mm ! From m/s to mm/s

        ! Calculate handling time, encounter rate, and probability of capture"""
        hand(j) = 0.264*10**(7.0151*(prey_LENGTH(j)/Larval_mm)) ! Walton 1992
        enc(j) = ((0.667*pi*(visual**3.)*travel + pi*(visual**2.)*sqrt(prey_LENGTH(j)**2.+ 2.*omega**2.)*travel*tau)*(MultiplyPrey*zoop(j))* ltr2mm3)
        ing(j) = (dt*enc(j)*pca(j)*prey_WGT(j)*ug2mg / (1 + hand(j)))*deltaH

       ! print *,"j:",j, " enc:",enc(j)*dt*deltaH, " mm:",Larval_mm, "pca:",pca(j), "SPpl:",prey_D(j)
        !print *,"j:",j, " ing:",ing(j), "suming:",sum(ing), "visual:",visual
    end do

    ! Calculate stomach fullness
    stomachFullness =  min(1.0,(Spre + sum(ing(:))/(Larval_wgt*gut_size)))
   ! write(*,*) "stomach fullness", stomachFullness
    suming = sum(ing(:))

    return

    end subroutine growth


double precision function n_dev (r_rand)
    Implicit None
    double precision U1, U2
    double precision r_rand, pi
    pi = 4.D0*DATAN(1.D0)
    r_rand = r_rand
    U1 = max (0.00001, r_rand)
    U2 = r_rand
    n_dev = Sqrt (-2.*dlog(U1)) * cos (2.*pi*U2)

    Return
  End Function n_dev

end module bioenergetics