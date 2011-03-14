Module zooplankton

  !
  ! ---------------------------------------------------------------
  ! This module calculates the module used to estimate the zooplankton
  ! abundance for each time step and depth level.
  !
  ! email : tk@trondkristiansen.com
  ! --------------------------------------------------------------------------------------------------------------------------------------------------
  !
  ! USAGE:
  !    var_array_rawXYZ[:,:,var_numberXYZ] = zooplankton.zooplankton.calculatezoo(np.asarray(accZoopArray, order='Fortran'),
  !                                                                                     NO3SM,NO3LG,
  !                                                                                      NH4SM,NH4LG,
  !                                                                                      TEMP,prey_D,
  !                                                                                      prey_WGT,numberOfDays,TT,DD,II)
  ! NOTE: that the dimensions of the arrays has to come last in the input order (TT,DD,II). This
  ! can be seen if you do: python; import zooplankton; print zooplankton.zooplankton.__doc__
  ! ---------------------------------------------------------------------------------------------------------------------------------------------------

  Implicit None

  ! f2py --verbose --fcompiler=intelem -c -m zooplankton zooplankton.f90 --f90flags="-no-heap-arrays"
  ! Also in .bash_profile it pays off substantially to only use :export FFLAGS="-m64 -O3 -fast -arch x86_64

Contains

    subroutine calculateZoo(accZoopArray,NO3SM,NO3LG,NH4SM,NH4LG,Temp,prey_D,prey_WGT,numberOfDays,TT,DD,II)
        implicit none

        Integer TT, DD, II
        Integer i, t, d, g

        double precision, dimension(TT,DD,II) :: accZoopArray
        double precision, dimension(TT,DD) :: NO3SM, NO3LG, NH4SM, NH4LG, Temp
        double precision, dimension(II) :: prey_D, prey_WGT

        double precision sizeMin, sizeMax, tm
        double precision, dimension(II):: prey_LENGTH, SDpl, SPpl, m
        double precision, dimension(II):: prey_AREA, prey_WIDTH

        double precision numberOfDays, Q10, Mo, Mqo, Mq, avgZ

        double precision, dimension(TT,DD,II) :: zoopArray
        double precision, dimension(TT,DD) :: NO3P_N, NO3P_C, NH4P_N, NH4P_C, NO3P, NH4P

!f2py intent(in,out,overwrite) accZoopArray
!f2py intent(in) NO3SM, NO3LG, NH4SM, NH4LG, Temp, prey_D, prey_WGT
!f2py intent(in) TT,DD,II, numberOfDays
!f2py intent(hide) PP

        !Convert production in N m-3 s-1 to mmoles C m-3 day-1
        !The production units are moles N m-3 sec-1, where N is nitrogen.  I'll
        !call this quantity PP_N.   So, if you wanted to convert to something
        !like mmoles C m-3 day-1 (PP_C), you would need to multiply PP_N by
        !a few conversions.
        !PP_C = PP_N * 6.625 moles C/mole N * 86400 sec/day * 1000 mmoles C/mole C
        !Note that 6.625 is 106/16 (i.e., Redfield C:N).

        NO3P_N = (NO3SM*0.25)**2 + NO3LG*0.25
        NH4P_N = (NH4SM*0.25)**2 + NH4LG*0.25

        !Convert production from nitrate to carbon
        NO3P_C = NO3P_N* 6.625 * 86400. * 1000.
        NH4P_C = NH4P_N* 6.625 * 86400. * 1000.

        !Convert from mmoles C m-3 d-1 to ug C m-3 d-1
        !(Carbon Atmoic Mass=12.0107)

        NO3P = NO3P_C * 12.0107 * 1000.
        NH4P = NH4P_C * 12.0107 * 1000.

        !Convert the production terms into biomass for each size class (SPpl).
        !The production terms for zooplankton is in ugCm-3-d-1 while growth
        !rate is d-1. The output is ug C m-3"""
        !New production

        do t=1,TT
            do d=1,DD
                do i=1,II
                    !Add new production
                    zoopArray(t,d,i) =(NO3P_C(t,d)) *numberOfDays
                    !Add old production
                    zoopArray(t,d,i) = zoopArray(t,d,i) + (NH4P(t,d) *numberOfDays)
                end do
            end do
        end do

        !Convert the units from ug C m-3 to ug dryweight m-3 assuming 50%
        !of zooplankton is carbon (need ref and perhaps other value...)
        zoopArray = zoopArray * 2.0

        ! STEP 2: Calculate the density function of prey of a given size according to Daewel et al. 2008"""

        sizeMin=100.        ! microm (/1000 to get mm)
        sizeMax=1600.

        !Calculate the sizes based on a range of 16 from 100 to 1600 micrograms
        do g=1,II
            prey_LENGTH(g)=sizeMin*(g)
        end do

        tm=0          !  total biomass
        SDpl(:)=0.0
        SPpl(:)=0.0
        m(:)=0.0
        prey_D(:)=0.0
        prey_AREA(:)=0.0

        prey_WGT = (/ 0.1, 0.25, 0.45, 1.0, 1.51, 2.76, 3.76, 6.0, 9.0, 13.24, 15.0, 17.0, 19.0, 23.13, 30.0, 45.0 /)
        prey_WIDTH= (/0.08, 0.1, 0.1, 0.1, 0.15, 0.17, 0.18, 0.2, 0.22, 0.25, 0.31, 0.35, 0.38, 0.38, 0.39, 0.40 /)

        do g=1,II
            SDpl(g) = 695.73 * exp(-0.0083*prey_LENGTH(g))
            m(g) = exp( 2.772 * log(prey_LENGTH(g)) - 7.476)
            tm = tm+m(g)
        end do

        !Here we calculate the realtive abundance (0-1) of total abundance for
        !each size group represented as a biomass percent"""
        do g=1,II
            SPpl(g) = (SDpl(g) * m(g))/tm
        end do

        prey_LENGTH= prey_LENGTH/1000. !micromm to mm

        !Calculate the image area (mm^2) of elongated prey
        do g=1,II
            prey_AREA(g)=0.75*(prey_LENGTH(g))*(prey_WIDTH(g))
        end do

        !Convert the biomass m-3 per stage to numbers per liter (for each time step)
        do t=1,TT
            do d=1,DD
                do i=1,II
                    zoopArray(t,d,i) = SPpl(i) * (zoopArray(t,d,i) /prey_WGT(i)) / 1000.

                    !Add the production to the total accumulated biomass"""
                    if (t .EQ. 1) then
                        accZoopArray(t,d,i)=zoopArray(t,d,i)
                    else
                        Q10=2.0
                        Mo=0.2
                        avgZ=50.

                        !Quadratic mortality term: Mq*Z_avg**2 = Mo*Z_avg"""
                        Mqo=Mo/(avgZ)

                        Mq=(Q10**((T-6.4)/10.) * Mqo)
                        Mq=Mq*accZoopArray(t-1,d,i)**2

                        accZoopArray(t,d,i)=exp(-Mq)*(accZoopArray(t-1,d,i)) + zoopArray(t,d,i)
                    end if
                    !print*,t,d,i,zoopArray(t,d,i),accZoopArray(t,d,i)
                end do
            end do
        end do
    end subroutine

end module zooplankton