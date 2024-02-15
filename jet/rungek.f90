       SUBROUTINE RUNGEK()
       use dimensiones
       use consadim
       use velocidades
       use variables
       use tiempo
       use right
       use jacobtools
       IMPLICIT NONE
       integer i,j,k,l,m

        dto=0.0
        dt=0.0001
        it=1
        DO WHILE (dto.le.tmax)
         CALL tstep()
         write(6,101)'it=',it,' ','dt=',dt,' ','dtotal=',dto
!______________________________________________________________________
!  1
         !write(6,*)'1'
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um0(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,m)=um0(i,j,k,m)+  &
                          dt*rs(i,j,k)  
            ENDDO
           ENDDO
          ENDDO
         ENDDO


         CALL varc_to_varnc()
         CALL frontera ()

         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo


!______________________________________________________________________
!  2
         write(6,*)'2'
         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,m)=0.75*um0(i,j,k,m)+  &
                          0.25*(um1(i,j,k,m)+ &
                          dt*rs(i,j,k))
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         

         CALL varc_to_varnc()
         CALL frontera ()
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo
         
!______________________________________________________________________
!  3
         write(6,*)'3'
         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,m)=(1./3.)*(um0(i,j,k,m)+   &
                                   2.*um1(i,j,k,m)+  &
                                   2.*dt*rs(i,j,k))
            ENDDO
           ENDDO
          ENDDO
         ENDDO


!         Write(6,*) 'FILTRADO'
         CALL ESTADISTICAS()
         CALL FILTRADO()
!         DO m=1,nd
!         CALL SOURCE(m)
!         ENDDO


         CALL varc_to_varnc()
         CALL frontera ()

!         do j=1,ny
!          write(6,*)um1(15,j,30,1),um1(15,j,30,1),temp(15,j,30)
!         enddo

         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um0(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo
         
         nver=mod(it/itr)
         if(nver.eq.0) then
          call MAXVAL()
         endif

         
         dto=dto+dt
         it=it+1
        ENDDO
        RETURN
 101    FORMAT(a3,i4,a1,a3,e12.4,a1,a7,e12.4)
        END SUBROUTINE RUNGEK

