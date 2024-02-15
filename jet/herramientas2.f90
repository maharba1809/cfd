!_______________________________________________________________________
      SUBROUTINE varc_to_varnc()
!_______________________________________________________________________
      use variables
      use velocidades
      use dimensiones
      use consadim
      IMPLICIT NONE
      integer :: i,j,k,l,m
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             u(i,j,k,1)=um1(i,j,k,2)/um1(i,j,k,1)
             u(i,j,k,2)=um1(i,j,k,3)/um1(i,j,k,1)
             u(i,j,k,3)=um1(i,j,k,4)/um1(i,j,k,1)

             um1(i,j,k,5)=(1.0/c2)*pres(i,j,k)+c1*um1(i,j,k,1)*  &
             (u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)

             temp(i,j,k)=c2*(um1(i,j,k,5)/um1(i,j,k,1)        &
                         -c1*(um1(i,j,k,2)**2+um1(i,j,k,3)**2 &
                         +um1(i,j,k,4)**2)/um1(i,j,k,1)**2.)
             pres(i,j,k)=um1(i,j,k,1)*temp(i,j,k)
            ENDDO
           ENDDO
          ENDDO

       RETURN
       END SUBROUTINE varc_to_varnc
!_______________________________________________________________________
      SUBROUTINE varnc_to_varc()
!_______________________________________________________________________
      use variables
      use velocidades
      use dimensiones
      use consadim
      IMPLICIT NONE
      integer :: i,j,k,l,m

       DO k=1,nz
        DO j=1,ny
         DO i=1,nx
          um0(i,j,k,1)=pres(i,j,k)/temp(i,j,k)
          um0(i,j,k,5)=(1.0/c2)*pres(i,j,k)+c1*um0(i,j,k,1)*  &
             (u(i,j,k,1)**2+u(i,j,k,2)**2+u(i,j,k,3)**2)
          um0(i,j,k,2)=u(i,j,k,1)*um0(i,j,k,1)
          um0(i,j,k,3)=u(i,j,k,2)*um0(i,j,k,1)
          um0(i,j,k,4)=u(i,j,k,3)*um0(i,j,k,1)
         ENDDO
        ENDDO
       ENDDO

       RETURN
       END SUBROUTINE varnc_to_varc



!-----------------------------------------------------------------------
      subroutine scatter(A,B,n)
!-----------------------------------------------------------------------

      use dimensiones

      implicit none

      integer :: n
      real, dimension(n) :: A,B

      B=A

      end subroutine scatter

!-----------------------------------------------------------------------
      subroutine scatter3d(A,B)
!-----------------------------------------------------------------------

      use dimensiones

      implicit none

      real, dimension(nx,ny,nz) :: A,B

      B=A

      end subroutine scatter3d

!-----------------------------------------------------------------------
      subroutine filtrado(itf)
!-----------------------------------------------------------------------

      use dimensiones
      use variables
      use consdernper_f
      use consderper_f
      use derivtools
      use jacobtools

      implicit none
      integer :: i,j,k,l,m,itf

       DO m=1,nd
        IF((itf.eq.1).or.(itf.eq.2)) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)
         ENDDO
        END DO
       END DO

        ELSEIF((itf.eq.3).or.(itf.eq.4)) THEN
        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)
          ENDDO
         END DO
        END DO

        ELSEIF((itf.eq.5).or.(itf.eq.6)) THEN
        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
           dw(k)=um1(i,j,k,m)*jbn(i,j,k,11)
          ENDDO
          call filtronper(nz,dw,aznpf,bznpf,cznpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,10)
          ENDDO
         END DO
        END DO
        ENDIF

        IF((itf.eq.3).or.(itf.eq.5)) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)
         ENDDO
        END DO
       END DO

        ELSEIF((itf.eq.1).or.(itf.eq.6)) THEN
        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)
          ENDDO
         END DO
        END DO

        ELSEIF((itf.eq.2).or.(itf.eq.4)) THEN
        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
          dw(k)=um1(i,j,k,m)*jbn(i,j,k,11)
         ENDDO
          call filtronper(nz,dw,aznpf,bznpf,cznpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,10)
          ENDDO
         END DO
        END DO
        ENDIF

        IF((itf.eq.6).or.(itf.eq.4)) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)
         ENDDO
        END DO
       END DO

        ELSEIF((itf.eq.2).or.(itf.eq.5)) THEN
        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)
          ENDDO
         END DO
        END DO

        ELSEIF((itf.eq.1).or.(itf.eq.3)) THEN
        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
           dw(k)=um1(i,j,k,m)*jbn(i,j,k,11)
          ENDDO
          call filtronper(nz,dw,aznpf,bznpf,cznpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,10)
          ENDDO
         END DO
        END DO
        ENDIF

       END DO

       RETURN
       END SUBROUTINE FILTRADO

!-----------------------------------------------------------------------
       SUBROUTINE FRONTERA ()
!-----------------------------------------------------------------------

       use dimensiones
       use variables
       use velocidades
       use consadim
       implicit none
      integer :: i,j,k,l,m
      real :: uu,vv,ww

       
       do j=1,ny
       do i=1,nx
        uu=um1(i,j,2,2)/um1(i,j,1,1)
        vv=um1(i,j,2,3)/um1(i,j,1,1)
        ww=0.0
        temp(i,j,1)=1.0
        pres(i,j,1)=pres(i,j,2)
        um1(i,j,1,1)=pres(i,j,1)/temp(i,j,1)
        um1(i,j,1,2)=uu*um1(i,j,2,1)
        um1(i,j,1,3)=vv*um1(i,j,2,1)
        um1(i,j,1,4)=ww*um1(i,j,2,1)
        u(i,j,1,1)=uu
        u(i,j,1,2)=vv
        u(i,j,1,3)=ww
        um1(i,j,1,5)=pres(i,j,1)/c2+   &
                     c1*um1(i,j,1,1)*(uu*uu+vv*vv+ww*ww)

        uu=um1(i,j,nz-1,2)/um1(i,j,nz-1,1)
        vv=um1(i,j,nz-1,3)/um1(i,j,nz-1,1)
        ww=0.0
        temp(i,j,nz)=1.2
        pres(i,j,nz)=pres(i,j,nz-1)
        um1(i,j,nz,1)=pres(i,j,nz)/temp(i,j,nz)
        um1(i,j,nz,2)=uu*um1(i,j,nz,1)
        um1(i,j,nz,3)=vv*um1(i,j,nz,1)
        um1(i,j,nz,4)=ww*um1(i,j,nz,1)
        u(i,j,nz,1)=uu
        u(i,j,nz,2)=vv
        u(i,j,nz,3)=ww
        um1(i,j,nz,5)=pres(i,j,nz)/c2+    &
                     c1*um1(i,j,nz,1)*(uu*uu+vv*vv+ww*ww)
       enddo
       enddo
       return
       end subroutine frontera

!_________________________________________________________________
      SUBROUTINE TP_MAKER()
!_________________________________________________________________
      use dimensiones
      use velocidades
      use mallagrid
      use vorticidad
      IMPLICIT NONE
      integer :: i,j,k,l,m
      !real, dimension (nx,ny,nz) :: buffer,wn
      real, allocatable,dimension (:,:,:) :: buffer,wn
      allocate(buffer(nx,ny,nz))
      allocate(wn(nx,ny,nz))
        DO k=1,nz
        DO j=1,ny
        DO i=1,nx
         wn(i,j,k)=(wx(i,j,k)**2.+wz(i,j,k)**2.+wy(i,j,k)**2.)**0.5
        ENDDO
        ENDDO
        ENDDO

!        OPEN(61,FILE='dataplt.dat',FORM='formatted')
!        write(61,*)'VARIABLES = "X", "Y", "Z", "U",  "T", "W"  &
!        , "V", "P","WN"'
!        write(61,*)'ZONE I=',NX,' J=',NY,' K=',NZ,' DATAPACKING=POINT'
!        DO k=1,nz
!        DO j=1,ny
!        DO i=1,nx
!        Write(61,*)x(i,j,k,1),x(i,j,k,2),x(i,j,k,3),u(i,j,k,1),  &
!           temp(i,j,k),u(i,j,k,3),u(i,j,k,2),pres(i,j,k),wn(i,j,k)
!        enddo
!        enddo
!        enddo
!        CLOSE (61)

      open(11,file='field_02.000',form='unformatted')
      write(11)nd
      write(11)Nx,Ny,Nz
      do k=1,nz
      do j=1,ny
      do i=1,nx
       buffer(i,j,k)=u(i,j,k,1)
      enddo
      enddo
      enddo
      write(11) buffer
      do k=1,nz
      do j=1,ny
      do i=1,nx
       buffer(i,j,k)=u(i,j,k,2)
      enddo
      enddo
      enddo
      write(11) buffer
      do k=1,nz
      do j=1,ny
      do i=1,nx
       buffer(i,j,k)=u(i,j,k,3)
      enddo
      enddo
      enddo
      write(11) buffer
      write(11) temp
      write(11) pres
      close(11)
      deallocate(buffer)
      deallocate(wn)
        RETURN
        end SUBROUTINE TP_MAKER
!_______________________________________________________
      SUBROUTINE MAXVAL()
!_______________________________________________________
   
      use dimensiones
      use velocidades
      use variables

      IMPLICIT NONE
      real umax,umin,zzz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer i,j,k,l,m

      umin=1000000.
      umax=0.00000001
      do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,1)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'umin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'umax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,2)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'vmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'vmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,3)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'wmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'wmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=pres(i,j,k)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'pmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'pmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=temp(i,j,k)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'tmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'tmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
      do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=um(i,j,k,1)
!if(j.eq.1.and.i.eq.1)write(*,*)k,zzz
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'rmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'rmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
 100  format(a4,a2,e12.4,a2,i3,a2,i3,a2,i3) 
      return
      end subroutine maxval
     

