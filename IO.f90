Module IO

Implicit None
  ! IO MODULE
  ! This module interpolates between time positions and depth levels to the
  ! current time and depth. So this is a 2D linear interpolation in
  ! time and space to get the excat value of several variables.
  !
  ! f2py --verbose --fcompiler=intelem -c -m IOdata io.f90 --f90flags="-no-heap-arrays"
  ! Also in .bash_profile it pays off substantially to only use :export FFLAGS="-m64 -O3 -fast -arch x86_64
  !
  ! Rewritten from Python module to Fortran module 29.01.2011
  ! Trond Kristiansen, tk@trondkristiansen.com
  ! www.trondkristiansen.com
  !
Contains

    subroutine getData(Tdata,Sdata,NH4SMdata,NH4LGdata,NO3SMdata,NO3LGdata,&
                    &CHLAdata,windX,windY,zoop,julian,julianIndex,julianFileA,&
                    &julianFileB,dz1,dz2,depthindex1,depthindex2,&
                    &inData,inDataXY,inDataZooplankton,event,&
                    &timeDim,depthDim,varDimXYZ,varDimXY,II)

        double precision Tdata,Sdata,NH4SMdata,NH4LGdata,NO3SMdata
        double precision TauX,TauY,NO3LGdata,CHLAdata
        double precision windX,windY,dz1,dz2
        double precision dwB, dwA
        integer julian,julianIndex,julianFileA, julianFileB
        integer depthindex1,depthindex2
        integer i, II
        integer timeDim, depthDim, varDimXYZ, varDimXY

        double precision, dimension(timeDim,depthDim,varDimXYZ):: inData
        double precision, dimension(timeDim,varDimXY):: inDataXY
        double precision, dimension(II):: zoop

        double precision, dimension(timeDim,depthDim,II):: inDataZooplankton
        character (len=32):: event

!f2py intent(in) julian,julianIndex,julianFileA,julianFileB,dz1,dz2,II,inDataZooplankton
!f2py intent(in) depthindex1,depthindex2,inData,inDataXY,event,timeDim,depthDim,varDim
!f2py intent(in,out,overwrite) Tdata,Sdata,NH4SMdata,NH4LGdata,NO3SMdata,NO3LGdata,CHLAdata,windX,windY,zoop

        !Calculate weights to use on input inData from file
        dwB = abs(julian) - abs(julianFileA)
        dwA = abs(julianFileB) - abs(julian)

        if (event=="REGULAR RUN" .OR. event=="TRISTAN RUN") then
            !Interpolate the values of temp, salt, u and v velocity in time to current julian date
            Tdata=((inData(julianIndex,depthindex1,1))*&
                &(dwA/(dwA+dwB))+(inData(julianindex+1,depthindex1,1))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,1))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,1))*&
                &(dwB/(dwA+dwB)))*dz2
            Sdata=((inData(julianIndex,depthindex1,2))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,2))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,2))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,2))*&
                &(dwB/(dwA+dwB)))*dz2
            TauX=((inDataXY(julianIndex,1))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,1))*&
                &(dwB/(dwA+dwB)))*dz1 +((inDataXY(julianIndex,1))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,1))*&
                &(dwB/(dwA+dwB)))*dz2
            TauY=((inDataXY(julianIndex,2))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,2))*&
                &(dwB/(dwA+dwB)))*dz1 +((inDataXY(julianIndex,2))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,2))*&
                &(dwB/(dwA+dwB)))*dz2
            !For regular runs, not all variables exists, but all will be
            !returned and therefore we have to define them as zero.
            NH4SMdata=0.0
            NH4LGdata=0.0
            NO3SMdata=0.0
            NO3LGdata=0.0
            CHLAdata=0.0

        else if (event=="ESM RUN") then
            !Interpolate the values of temp, salt, u and v velocity in time to current julian date.
            !Note that the indices for 2D and 3D values are different and must be set accordingly to
            !the order the 2D and 3D variables are stored into arrays. The order is defined by the order they appear
            !in vars list defined in init funtion.
            Tdata=((inData(julianIndex,depthindex1,1))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,1))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,1))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,1))*&
                &(dwB/(dwA+dwB)))*dz2
            Sdata=((inData(julianIndex,depthindex1,2))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,2))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,2))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,2))*&
                &(dwB/(dwA+dwB)))*dz2
            NH4SMdata=((inData(julianIndex,depthindex1,3))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,3))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,3))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,3))*&
                &(dwB/(dwA+dwB)))*dz2
            NH4LGdata=((inData(julianIndex,depthindex1,4))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,4))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,4))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,4))*&
                &(dwB/(dwA+dwB)))*dz2
            NO3SMdata=((inData(julianIndex,depthindex1,5))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,5))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,5))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,5))*&
                &(dwB/(dwA+dwB)))*dz2
            NO3LGdata=((inData(julianIndex,depthindex1,6))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,6))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,6))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,6))*&
                &(dwB/(dwA+dwB)))*dz2
            CHLAdata=((inData(julianIndex,depthindex1,7))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex1,7))*&
                &(dwB/(dwA+dwB)))*dz1 +((inData(julianIndex,depthindex2,7))*&
                &(dwA/(dwA+dwB))+(inData(julianIndex+1,depthindex2,7))*&
                &(dwB/(dwA+dwB)))*dz2
            TauX=((inDataXY(julianIndex,1))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,1))*&
                &(dwB/(dwA+dwB)))*dz1 +((inDataXY(julianIndex,1))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,1))*&
                &(dwB/(dwA+dwB)))*dz2
            TauY=((inDataXY(julianIndex,2))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,2))*&
                &(dwB/(dwA+dwB)))*dz1 +((inDataXY(julianIndex,2))*&
                &(dwA/(dwA+dwB))+(inDataXY(julianIndex+1,2))*&
                &(dwB/(dwA+dwB)))*dz2
            do i=1,II
                zoop(i)=((inDataZooplankton(julianIndex,depthindex1,i))*&
                &(dwA/(dwA+dwB))+(inDataZooplankton(julianIndex+1,depthindex1,i))*&
                &(dwB/(dwA+dwB)))*dz1 +((inDataZooplankton(julianIndex,depthindex2,i))*&
                &(dwA/(dwA+dwB))+(inDataZooplankton(julianIndex+1,depthindex2,i))*&
                &(dwB/(dwA+dwB)))*dz2
                !print*,"1",i,inDataZooplankton(julianIndex,depthindex1,i),inDataZooplankton(julianIndex+1,depthindex1,i)
                !print*,"2",i,inDataZooplankton(julianIndex,depthindex2,i),inDataZooplankton(julianIndex+1,depthindex2,i)
            end do
        else
            print*,"No specific type of run defined in getData (IO.f90)"




        call convertStressToWind(TauX,TauY,windX,windY)

        end if

        return

    end subroutine getData

    subroutine convertStressToWind(TauX,TauY,windX,windY)
        !Convert surface stress to wind velocity for use in turbulence calculations of
        !encounter rate between larvae and prey. We assume a drag coefficient
        !with the shape Cd=0.44+0.063U (where U is wind at 10 m).
        !In SODA TauX and TauY are in N/m**2. This is equivavlent to
        !N=kgm/s**2 and therefore kg/ms**2. The unit for air density is kg/m**3 so
        !after doing the quadratic estimate of U we get wind velocity with units m/s.

        !Update: changed drag coefficient top constant = 1.2e-3 and wind density to 1.3. This sfunction also
        !disregards direction of wind as we are only concerned with the scalar value."""

        double precision TauX, TauY, windX, windY
        windX=abs(sqrt(abs(TauX)/(1.3*1.2e-3)))
        windY=abs(sqrt(abs(TauY)/(1.3*1.2e-3)))

        return
    end subroutine convertStressToWind

end module IO