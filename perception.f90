Module perception

  !
  ! ---------------------------------------------------------------
  ! This module calculates the preceptive distance a cod larva
  ! of a certain weight and size is able to detect prey.
  ! Modified from Fiksen & MacKenzie (2002).

  ! Trond Kristiansen, 27.03.2004
  ! Last revision: 25.04.2004, 25.03.2005 Fixed a bug in surface_light routine
  ! Now the routine uses the correct number of days since January first YYYY
  ! to calculate the surface light(day_of_year,hour,latitude). Also changed the
  ! parameters to be correct according to Aksnes & Utne 1997
  !
  ! Converted into a python module using f2py on January 25th 2010 on train
  ! from Oslo to Skien.
  ! email : Trond.Kristiansen@imr.no
  ! ---------------------------------------------------------------

  Implicit None
  !

  ! f2py --verbose --fcompiler=intelem -c -m perception perception.f90 --f90flags="-no-heap-arrays"
Contains
  !
  Subroutine getr (r, c, C0, Ap, Vc, KE, EB, IER)
    ! USAGE:  r, IER = perception.perception.getr(r,beamAttCoeff,contrast,calanus_Area,Em,Ke_larvae,Eb, IER)
    ! where r is the visual range of prey item (mm),
    ! the beam attenuation coefficient (c) is usually 3*attenuation coefficient (~0.4 m^-1)
    ! (note if input units are in mm you have to divide attenuation coefficient with 1000 to
    ! get per mm instead of per m, ~0.4 -> 0.0004)
    ! contrast (C0) is the contrast of the prey to the background (~0.4),calanus_Area
    ! (Ap) is the size of the prey (the image area = 0.75*length*width in mm^2),
    ! Em is visual sensitivity of the predator
    ! (Em = Larval_mm**2/(contrast*0.1*0.2*0.75) - Size-specific sensitivity of the visual system (Fiksen,02))
    ! for detecting a pre item of width 0.1 and length 0.2 mm with shape elongated (0.75)
    ! Input of larval length (larval_mm) is in units mm. Ke_larvae is the visual threshold for detecting prey
    ! at low light levels (umol/m2/s-1) which is usually set to Ke_larvae=1. Eb is the light level at the
    ! current depth level (use calclight.f90 to estimate light).
    ! Edited by Trond Kristiansen, 26.01.2010.
    !
    !     -------------------------------------------------------------
    !     Obtain visual range by solving the non-linear equation
    !     by means of Newton-Raphson iteration and derivation in
    !     subroutine DERIV. Initial value is calculated in EASYR.
    !	  The calculation is based on the model described in Aksnes &
    !	  Utne (1997) Sarsia 83:137-147.
    !
    !	  Programmed and tested 29.01.01 Dag L Aksnes

    Implicit None

    Real r, c, C0, Vc, Ap, KE, EB
    Real EPS, RST, TOL, TOLF, F1, FDER, DX
    Integer I, IEND, AS, IER
    !
    !     Input parameters
    !		RST		: start value of r calculated by EASYR
    !		c		: beam attenuation koefficient (m-1)
    !		C0		: prey inherent contrast
    !		Ap		: prey area (m^2)
    !		Vc		: parameter characterising visual capacity (d.l.)
    !				  this parameter is denoted E' in Aksnes & Utne
    !	   	Ke	    : saturation parameter (uE m-2 s-1)
    !	    Eb	    : background irradiance at depth DEPTH
    !
    !    Output parameters
    !        F1     : function value of equation in DERIV
    !        FDER   : the derivative of the function
    !        r      : the predator's visual range (when F1=0)
    !	   IER	  : = 1, No convergence after IEND steps. Error return.
    !				= 2, Return in case of zero divisor.
    !				= 3, r out of allowed range (negative)
    !				= 0, valid r returned
    !
    !

!f2py intent(in,overwrite) c, C0, Ap, Vc, KE, EB, IER
!f2py intent(in,out,overwrite) r

    !
    !
    !     Initial guess of visual range (RST)
    Call EASYR (RST, C0, Ap, Vc, KE, EB)
    !
    !
    !     Upper boundary of allowed error of visual range
    EPS = .0001
    !     Maximum number of iteration steps
    IEND = 200
    !
    !     Prepare iteration
    r = RST
    TOL = r
    Call DERIV (r, F1, FDER, c, C0, Ap, Vc, KE, EB)
    TOLF = 100. * EPS
    !
    !     Start iteration loop
    Do I = 1, IEND
       If (F1) 1, 7, 1
       !
       !        Equation is not satisfied by r
1      If (FDER) 2, 8, 2
       !
       !        Iteration is possible
2      DX = F1 / FDER
       r = r - DX
       !
       !        Test on allowed range
       If (r .Lt. 0.) Go To 9
       !
       TOL = r
       Call DERIV (r, F1, FDER, c, C0, Ap, Vc, KE, EB)
       !
       !        Test on satisfactory accuracy
       TOL = EPS
       AS = Abs (r)
       If (AS-1.) 4, 4, 3
3      TOL = TOL * AS
4      If (Abs(DX)-TOL) 5, 5, 6
5      If (Abs(F1)-TOLF) 7, 7, 6
6      Continue
       !
       !     No convergence after IEND steps. Error return.
       IER = 1
7      Return
       !     Return in case of zero divisor
8      IER = 2
       Return
       !     r out of allowed range (negative)
9      IER = 3
       Return
       !
    End Do
  End Subroutine getr
  !
  !     -------------------------------------------------------------------
  Subroutine easyr (r, C0, Ap, Vc, KE, EB)
    !     -------------------------------------------------------------------
    !     Obtain a first estimate of visual range by using a simplified
    !     expression of visual range
    Implicit None
    Real r, C0, Ap, Vc, KE, EB
    Real R2
    !
    !     Se calling routine for explanation of parameters
    !
    R2 = Abs (C0) * Ap * Vc * (EB/(KE+EB))
    r = Sqrt (R2)
    !write(*,*), "easyR",r
    Return
  End Subroutine easyr
  !
  !
  !     -------------------------------------------------------
  Subroutine deriv (r, F1, FDER, c, C0, Ap, Vc, KE, EB)
    !     -------------------------------------------------------

    !     Derivation of equation for visual range of a predator
    Implicit None
    Real r, c, C0, Ap, Vc, KE, EB
    Real FR1, FR2, F1, FDER
    !
    !     Input parameters
    !        Se explanation in calling routine
    !
    !    Output parameters
    !        F1     : function value of equation in DERIV
    !        FDER   : the derivative of the function
    !        r      : the predator's visual range (when F1=0)
    !
    !	   The function and the derivative is calculated on the basis of the
    !	   log-transformed expression
    !
    !
    FR2 = Log (Abs(C0)*Ap*Vc)
    FR1 = Log (((KE+EB)/EB)*r*r*Exp(c*r))
    F1 = FR1 - FR2
    FDER = c + 2. / r
    !
    Return
  End Subroutine deriv

End Module perception
