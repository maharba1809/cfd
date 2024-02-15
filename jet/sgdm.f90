!_________________________________________________________________
       SUBROUTINE INISGDM()
!_________________________________________________________________
      use sgdmodel
       use dimensiones
       use mallagrid
       IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: dxmax,dymax,dzmax

       isgm=2
       if(isgm.eq.1) then
        amu0=0.063
       elseif(isgm.eq.2) then
        amu0=0.1781
        cutoff=COS(20.0*ATAN(1.0)/45.0) 
       endif
       do k=1,nz
        do j=1,ny
         do i=2,nx-1
          dxmsgd(i,j,k)=x(i,j,k,1)-x(i-1,j,k,1)   
          dxpsgd(i,j,k)=x(i+1,j,k,1)-x(i,j,k,1)
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do j=1,ny
         dxpsgd(1,j,k)=(x(2,j,k,1)-x(1,j,k,1))
         dxmsgd(1,j,k)=(x(2,j,k,1)-x(1,j,k,1))
         dxmsgd(nx,j,k)=(x(nx,j,k,1)-x(nx-1,j,k,1))
         dxpsgd(nx,j,k)=(x(nx,j,k,1)-x(nx-1,j,k,1))
        ENDDO
       ENDDO

       do k=1,nz
        do j=2,ny-1
         do i=1,nx
          dymsgd(i,j,k)=x(i,j,k,2)-x(i,j-1,k,2)
          dypsgd(i,j,k)=x(i,j+1,k,2)-x(i,j,k,2)
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do i=1,nx
         dypsgd(i,1,k)=(x(i,2,k,2)-x(i,1,k,2))
         dymsgd(i,1,k)=(x(i,2,k,2)-x(i,1,k,2))
         dymsgd(i,ny,k)=(x(i,ny,k,2)-x(i,ny-1,k,2))
         dypsgd(i,ny,k)=(x(i,ny,k,2)-x(i,ny-1,k,2))
        ENDDO
       ENDDO

       do k=2,nz-1
        do j=1,ny
         do i=1,nx
          dzmsgd(i,j,k)=x(i,j,k,3)-x(i,j,k-1,3)
          dzpsgd(i,j,k)=x(i,j,k+1,3)-x(i,j,k,3)
         ENDDO
        ENDDO
       ENDDO
       do j=1,ny
        do i=1,nx
         dzpsgd(i,j,1)=(x(i,j,2,3)-x(i,j,1,3))
         dzmsgd(i,j,1)=(x(i,j,2,3)-x(i,j,1,3))
         dzmsgd(i,j,nz)=(x(i,j,nz,3)-x(i,j,nz-1,3))
         dzpsgd(i,j,nz)=(x(i,j,nz,3)-x(i,j,nz-1,3))
        ENDDO
       ENDDO
       do K=1,nz
        do j=1,ny
         do i=1,nx
         dzmax=max(dzmsgd(i,j,k),dzpsgd(i,j,k))
         dymax=max(dymsgd(i,j,k),dypsgd(i,j,k))
         dxmax=max(dxmsgd(i,j,k),dxpsgd(i,j,k))
         dxyzsgd(i,j,k)=(dxmax*dymax*dzmax)**(1./3.)
         enddo
        enddo
       enddo
       do K=1,nz
        do j=1,ny
         do i=1,nx
         dxpsgd(i,j,k)=(dxyzsgd(i,j,k)/dxpsgd(i,j,k))**(1.0/3.0)
         dxmsgd(i,j,k)=(dxyzsgd(i,j,k)/dxmsgd(i,j,k))**(1.0/3.0)
         dypsgd(i,j,k)=(dxyzsgd(i,j,k)/dypsgd(i,j,k))**(1.0/3.0)
         dymsgd(i,j,k)=(dxyzsgd(i,j,k)/dymsgd(i,j,k))**(1.0/3.0)
         dzpsgd(i,j,k)=(dxyzsgd(i,j,k)/dzpsgd(i,j,k))**(1.0/3.0)
         dymsgd(i,j,k)=(dxyzsgd(i,j,k)/dzmsgd(i,j,k))**(1.0/3.0)
         enddo
        enddo
       enddo




       RETURN
       END SUBROUTINE INISGDM

!_________________________________________________________________
       SUBROUTINE viscosidad_t ()
!_________________________________________________________________

       use dimensiones
       use sgdmodel
       use variables
       IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: umax,zzz
      INTEGER :: kmax,jmax,imax

      CALL f_structura ()

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         amut(i,j,k)=amu0*SQRT(amut(i,j,k))*dxyzsgd(i,j,k)*um(i,j,k,1)
        END DO
       END DO
      END DO

      if(isgm.eq.2) then
       CALL vort ()
       CALL selectivo ()
      endif

      umax=-1.e10
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
          zzz=amut(i,j,k)
          IF (zzz.gt.umax) THEN
           umax=zzz
           imax=i
           jmax=j
           kmax=k
          ENDIF
         END DO
        END DO
      END DO

      !PRINT *,'nutmax,i,j,k= ',umax,imax,jmax,kmax


      RETURN
      end SUBROUTINE VISCOSIDAD_t


!_________________________________________________________________
      SUBROUTINE f_structura ()
!_________________________________________________________________

       use dimensiones
       use sgdmodel
       use variables
       use velocidades
       IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: u1,v1,w1,u2,v2,w2,u3,v3,w3
      real :: u4,v4,w4,u5,v5,w5,u6,v6,w6

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx-1
         u2=(u(i+1,j,k,1)-u(i,j,k,1))*dxpsgd(i,j,k)
         v2=(u(i+1,j,k,2)-u(i,j,k,2))*dxpsgd(i,j,k)
         w2=(u(i+1,j,k,3)-u(i,j,k,3))*dxpsgd(i,j,k)
         amut(i,j,k)=u2*u2+v2*v2+w2*w2
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
         u2=(u(nx,j,k,1)-u(nx-1,j,k,1))*dxpsgd(nx,j,k)
         v2=(u(nx,j,k,2)-u(nx-1,j,k,2))*dxpsgd(nx,j,k)
         w2=(u(nx,j,k,3)-u(nx-1,j,k,3))*dxpsgd(nx,j,k)
         amut(nx,j,k)=u2*u2+v2*v2+w2*w2
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        DO i=2,nx
         u1=(u(i-1,j,k,1)-u(i,j,k,1))*dxmsgd(i,j,k)
         v1=(u(i-1,j,k,2)-u(i,j,k,2))*dxmsgd(i,j,k)
         w1=(u(i-1,j,k,3)-u(i,j,k,3))*dxmsgd(i,j,k)
         amut(i,j,k)=amut(i,j,k)+u1*u1+v1*v1+w1*w1
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        u1=(u(2,j,k,1)-u(1,j,k,1))*dxmsgd(1,j,k)
        v1=(u(2,j,k,2)-u(1,j,k,2))*dxmsgd(1,j,k)
        w1=(u(2,j,k,3)-u(1,j,k,3))*dxmsgd(1,j,k)
        amut(1,j,k)=amut(1,j,k)+u1*u1+v1*v1+w1*w1
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny-1
        DO i=1,nx
         u3=(u(i,j+1,k,1)-u(i,j,k,1))*dypsgd(i,j,k)
         v3=(u(i,j+1,k,2)-u(i,j,k,2))*dypsgd(i,j,k)
         w3=(u(i,j+1,k,3)-u(i,j,k,3))*dypsgd(i,j,k)
         amut(i,j,k)=amut(i,j,k)+u3*u3+v3*v3+w3*w3
        END DO
       END DO
      END DO
      DO k=1,nz
       DO i=1,nx
        u3=(u(i,ny,k,1)-u(i,ny-1,k,1))*dypsgd(i,ny,k)
        v3=(u(i,ny,k,2)-u(i,ny-1,k,2))*dypsgd(i,ny,k)
        w3=(u(i,ny,k,3)-u(i,ny-1,k,3))*dypsgd(i,ny,k)
        amut(i,ny,k)=amut(i,ny,k)+u3*u3+v3*v3+w3*w3
       END DO
      END DO
      DO k=1,nz
       DO j=2,ny
        DO i=1,nx
         u4=(u(i,j-1,k,1)-u(i,j,k,1))*dymsgd(i,j,k)
         v4=(u(i,j-1,k,2)-u(i,j,k,2))*dymsgd(i,j,k)
         w4=(u(i,j-1,k,3)-u(i,j,k,3))*dymsgd(i,j,k)
         amut(i,j,k)=amut(i,j,k)+u4*u4+v4*v4+w4*w4
        END DO
       END DO
      END DO
      DO k=1,nz
       DO i=1,nx
        u4=(u(i,1,k,1)-u(i,2,k,1))*dymsgd(i,1,k)
        v4=(u(i,1,k,2)-u(i,2,k,2))*dymsgd(i,1,k)
        w4=(u(i,1,k,3)-u(i,2,k,3))*dymsgd(i,1,k)
        amut(i,1,k)=amut(i,1,k)+u4*u4+v4*v4+w4*w4
       END DO
      END DO

      DO k=1,nz-1
       DO j=1,ny
        DO i=1,nx
         u5=(u(i,j,k+1,1)-u(i,j,k,1))*dzpsgd(i,j,k)
         v5=(u(i,j,k+1,2)-u(i,j,k,2))*dzpsgd(i,j,k)
         w5=(u(i,j,k+1,3)-u(i,j,k,3))*dzpsgd(i,j,k)
         amut(i,j,k)=amut(i,j,k)+u5*u5+v5*v5+w5*w5
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        u5=(u(i,j,nz,1)-u(i,j,nz-1,1))*dzpsgd(i,j,nz)
        v5=(u(i,j,nz,2)-u(i,j,nz-1,2))*dzpsgd(i,j,nz)
        w5=(u(i,j,nz,3)-u(i,j,nz-1,3))*dzpsgd(i,j,nz)
        amut(i,j,nz)=amut(i,j,nz)+u5*u5+v5*v5+w5*w5
       END DO
      END DO
      DO k=2,nz
       DO j=1,ny
        DO i=1,nx
         u6=(u(i,j,k-1,1)-u(i,j,k,1))*dzmsgd(i,j,k)
         v6=(u(i,j,k-1,2)-u(i,j,k,2))*dzmsgd(i,j,k)
         w6=(u(i,j,k-1,3)-u(i,j,k,3))*dzmsgd(i,j,k)
         amut(i,j,k)=amut(i,j,k)+u6*u6+v6*v6+w6*w6
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        u6=(u(i,j,1,1)-u(i,j,2,1))*dzmsgd(i,j,1)
        v6=(u(i,j,1,2)-u(i,j,2,2))*dzmsgd(i,j,1)
        w6=(u(i,j,1,3)-u(i,j,2,3))*dzmsgd(i,j,1)
        amut(i,j,1)=amut(i,j,1)+u6*u6+v6*v6+w6*w6
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         amut(i,j,k)=amut(i,j,k)/6.e0
        END DO
       END DO
      END DO

      RETURN
      END SUBROUTINE f_structura


!_________________________________________________________________
      SUBROUTINE VORT ()
!_________________________________________________________________

       use dimensiones
       use velocidades
       use vorticidad
       use jacobtools
       use derivvel
       IMPLICIT NONE
       integer :: i,j,k,l,m

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         wx(i,j,k)=jbn(i,j,k,5)*dvel(i,j,k,8)-  &
                   jbn(i,j,k,9)*dvel(i,j,k,6)
         wy(i,j,k)=jbn(i,j,k,9)*dvel(i,j,k,3)-  &
                   jbn(i,j,k,1)*dvel(i,j,k,7)
         wz(i,j,k)=jbn(i,j,k,1)*dvel(i,j,k,4)-  &
                   jbn(i,j,k,5)*dvel(i,j,k,2)  
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE VORT

!_________________________________________________________________
       SUBROUTINE SELECTIVO ()
!_________________________________________________________________

       use dimensiones
       use sgdmodel
       use variables
       use vorticidad
       IMPLICIT NONE
      integer :: i,j,k,l,m,jm,jp,im,ip,kp,km
      real :: wxn,wyn,wzn,ps1,ps2,ps3,wn,s1




      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         ip=i+1
         im=i-1
         jp=j+1
         jm=j-1
         kp=k+1
         km=k-1
         IF (i.eq.nx) ip=nx
         IF (i.eq.1)  im=1
         IF (j.eq.ny) jp=ny
         IF (j.eq.1)  jm=1
         IF (k.eq.nz) kp=nz
         IF (k.eq.1)  km=1
         wxn=wx(ip,j,k)+wx(im,j,k)+wx(i,j,kp)+wx(i,j,km) &
           +wx(i,jp,k)+wx(i,jm,k)
         wyn=wy(ip,j,k)+wy(im,j,k)+wy(i,j,kp)+wy(i,j,km) &
           +wy(i,jp,k)+wy(i,jm,k)
         wzn=wz(ip,j,k)+wz(im,j,k)+wz(i,j,kp)+wz(i,j,km) &
           +wz(i,jp,k)+wz(i,jm,k)
         WN=wx(i,j,k)**2+wy(i,j,k)**2+wz(i,j,k)**2
         WN=SQRT(WN*(wxn**2+wyn**2+wzn**2))
         PS1=wx(i,j,k)*wxn
         PS2=wy(i,j,k)*wyn
         PS3=wz(i,j,k)*wzn
         S1=(PS1+PS2+PS3)/WN
         IF (S1.gt.cutoff) THEN
           amut(i,j,k)=0.0
         ENDIF
        END DO
       END DO
      END DO

       RETURN
       END SUBROUTINE SELECTIVO
