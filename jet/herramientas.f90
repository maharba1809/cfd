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
             temp(i,j,k)=c2*(um1(i,j,k,5)/um1(i,j,k,1)        &
                         -c1*(um1(i,j,k,2)**2+um1(i,j,k,3)**2 &
                         +um1(i,j,k,4)**2)/um1(i,j,k,1)**2.)
             pres(i,j,k)=um1(i,j,k,1)*temp(i,j,k)
             u(i,j,k,1)=um1(i,j,k,2)/um1(i,j,k,1)
             u(i,j,k,2)=um1(i,j,k,3)/um1(i,j,k,1)
             u(i,j,k,3)=um1(i,j,k,4)/um1(i,j,k,1)
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
      subroutine filtrado()
!-----------------------------------------------------------------------

      use dimensiones
      use variables
      use consdernper_f
      use consderper_f
      use derivtools
      use jacobtools

      implicit none
      integer :: i,j,k,l,m

       DO m=1,nd
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
         call filtroper(nx,du,axpf,bxpf,cxpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)*jbn(i,j,k,11)
         ENDDO
        END DO
       END DO

        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)
          ENDDO
          call filtroper(ny,dv,aypf,bypf,cypf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)
          ENDDO
         END DO
        END DO

        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
           dw(k)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(nz,dw,aznpf,bznpf,cznpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

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

       
       do j=1,ny
       do i=1,nx
        um1(i,j,1,2)=0.0
        um1(i,j,nz,2)=0.0
        um1(i,j,1,3)=0.0
        um1(i,j,nz,3)=0.0
        um1(i,j,1,4)=0.0
        um1(i,j,nz,4)=0.0
        u(i,j,nz,1)=0.0
        u(i,j,1,1)=0.0
        u(i,j,nz,2)=0.0
        u(i,j,1,2)=0.0
        u(i,j,nz,3)=0.0
        u(i,j,1,3)=0.0

        temp(i,j,1)=1.0
        temp(i,j,nz)=1.0
        pres(i,j,1)=pres(i,j,2)
        pres(i,j,nz)=pres(i,j,nz-1)
        um1(i,j,1,1)=pres(i,j,1)/temp(i,j,1)
        um1(i,j,nz,1)=pres(i,j,nz)/temp(i,j,nz)
        um1(i,j,1,5)=pres(i,j,1)/c2
        um1(i,j,nz,5)=pres(i,j,nz)/c2
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
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real, dimension (nx,ny,nz) :: buffer,wn

        DO k=1,nz
        DO j=1,ny
        DO i=1,nx
         wn(i,j,k)=(wx(i,j,k)**2.+wy(i,j,k)**2.+wz(i,j,k)**2.)**0.5
        ENDDO
        ENDDO
        ENDDO

        OPEN(61,FILE='dataplt.dat',FORM='formatted')
        write(61,*)'VARIABLES = "X", "Y", "Z", "U",  "T", "W"  &
        , "V", "P","WN"'
        write(61,*)'ZONE I=',NX,' J=',NY,' K=',NZ,' DATAPACKING=POINT'
        DO k=1,nz
        DO j=1,ny
        DO i=1,nx
         Write(61,*)x(i,j,k,1),x(i,j,k,2),x(i,j,k,3),u(i,j,k,1),  &
           temp(i,j,k),u(i,j,k,3),u(i,j,k,2),pres(i,j,k),wn(i,j,k)
        enddo
        enddo
        enddo
        CLOSE (61)

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
        RETURN
        end SUBROUTINE TP_MAKER


