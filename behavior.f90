Module behavior

Implicit None
  ! BEHAVIOR MODULE
  !
  ! f2py --verbose --fcompiler=intelem -c -m behavior behavior.f90 --f90flags="-no-heap-arrays"
  ! Also in .bash_profile it pays off substantially to only use :export FFLAGS="-m64 -O3 -fast -arch x86_64
  !
  ! Rewritten from Python module to Fortran module 29.01.2011, 07.11.2011
  ! Trond Kristiansen, tk@trondkristiansen.com
  ! www.trondkristiansen.com
  !
Contains

    subroutine getBehavior(stomachFullness,F,mortality,length,depth,optDepth,oldFitness)

        double precision stomachFullness, F, mortality, length
        double precision optDepth, oldFitness, depth
        double precision beta, T, vector, fitness

!f2py intent (in) stomachFullness,F,mortality,length,depth
!f2py intent (in,out,overwrite) optDepth,oldFitness

        !Rule 4 of Behavioral Ecology paper - Kristiansen et al. 2009"""
        T=min(0.9,0.3+1000.0*(1+length*exp(length))**(-1))

        if (stomachFullness .GT. T) then
            beta   = 7.0
            vector = ((stomachFullness-T)/(1.0-T))**beta
        else
            vector=0.0
        end if

        fitness=(((1-vector)*F) - vector*mortality)

        if (fitness .GT. oldFitness) then
         !   print "Change: Vector %s new depth %s vs old depth %s : new fit %s and old %s "%(vector,optDepth,depth,fitness,oldFitness)
            optDepth=depth
            oldFitness=fitness
        end if

        return

    end subroutine getBehavior

end module behavior