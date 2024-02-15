
       SUBROUTINE INITSTEP()
       use tiempo
       use dimensiones
       use mallagrid
       IMPLICIT NONE
      integer :: i,j,k,l,m

       up=0.5
       down=1.2
       do k=1,nz
        do j=1,ny
         do i=2,nx-1
          dx(i,j,k)=MIN(x(i,j,k,1)-x(i-1,j,k,1),   &
                        x(i+1,j,k,1)-x(i,j,k,1))
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do j=1,ny
         dx(1,j,k)=(x(2,j,k,1)-x(1,j,k,1))
         dx(nx,j,k)=(x(nx,j,k,1)-x(nx-1,j,k,1))
        ENDDO
       ENDDO

       do k=1,nz
        do j=2,ny-1
         do i=1,nx
          dy(i,j,k)=MIN(x(i,j,k,2)-x(i,j-1,k,2),   &
                        x(i,j+1,k,2)-x(i,j,k,2))
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do i=1,nx
         dy(i,1,k)=(x(i,2,k,2)-x(i,1,k,2))
         dy(i,ny,k)=(x(i,ny,k,2)-x(i,ny-1,k,2))
        ENDDO
       ENDDO

       do k=2,nz-1
        do j=1,ny
         do i=1,nx
          dz(i,j,k)=MIN(x(i,j,k,3)-x(i,j,k-1,3),   &
                        x(i,j,k+1,3)-x(i,j,k,3))
         ENDDO
        ENDDO
       ENDDO
       do j=1,ny
        do i=1,nx
         dz(i,j,1)=(x(i,j,2,3)-x(i,j,1,3))
         dz(i,j,nz)=(x(i,j,nz,3)-x(i,j,nz-1,3))
        ENDDO
       ENDDO


       do k=1,nz
        do j=1,ny
         do i=1,nx
          dy(i,j,k)=1./dy(i,j,k)
          dx(i,j,k)=1./dx(i,j,k)
          dz(i,j,k)=1./dz(i,j,k)
         ENDDO
        ENDDO
       ENDDO

       RETURN
       END SUBROUTINE INITSTEP
     

       SUBROUTINE TSTEP()
       use tiempo
       use dimensiones
       use velocidades
       use consadim
       use flow
       IMPLICIT NONE
      integer :: i,j,k,l,m
       real :: c


       olddt = dt
       cfla=0.
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
          c=SQRT(gamma*c3*temp(i,j,k))
          cfla = MAX(cfla,dx(i,j,k)*(ABS(u(i,j,k,1))+c), &
                          dy(i,j,k)*(ABS(u(i,j,k,2))+c), &
                          dz(i,j,k)*(ABS(u(i,j,k,3))+c))
        END DO
       END DO
      END DO
      cflm = cfla
      cflm = olddt*cflm
!
!       ----------------------------------------------
!       RNUEVO PASO DE TIEMPO
!       ----------------------------------------------
!
      dt = olddt/cflm*cfl
      diffdt = dt - olddt
      dt = olddt + 0.5*(up+down)*diffdt &
                + 0.5*(up-down)*ABS(diffdt)
       RETURN
       END SUBROUTINE TSTEP

