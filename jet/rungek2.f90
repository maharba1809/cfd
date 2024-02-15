       SUBROUTINE RUNGEK()
       use dimensiones
       use consadim
       use velocidades
       use variables
       use tiempo
       use right
       use jacobtools
       use fieldrec
       use vieu
       use fronterapl
       IMPLICIT NONE
       integer i,j,k,l,m,agrava,avieu
       real xm,pe,re

        dto=0.0
        dt=0.0001
        it=1
        itf=1
        igrava=0
        ivieu=1
        dgrava=dgrava1
        tvieu=tvieu1
        i_sample=0
!        CALL MAXVAL()
        CALL GRAVA ()
        DO WHILE (dto.le.tmax)
         CALL tstep()
     write(6,*)'it=',it,' ','dt=',dt,' ','dtotal=',dto,'dgrava',dgrava
!     write(6,101)'it=',it,' ','dt=',dt,' ','dtotal=',dto,'dgrava',dgrava
!______________________________________________________________________
!  1
         write(6,*)'1'
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um0(i,j,k,m)
         enddo
         enddo
         enddo
         enddo
        CALL MAXVAL()

         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         CALL esponja()
         CALL frontera_x()
         CALL frontera_z()
         CALL frontera_y()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx

             um1(i,j,k,m)=um0(i,j,k,m)+  &
                          dt*(rs(i,j,k)+rsv(i,j,k))
            ENDDO
           ENDDO
          ENDDO
          CALL FXSOURCE (1,m)
          CALL FZSOURCE (1,m)
          CALL FYSOURCE (1,m)
         ENDDO

          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
! presion pe
             pe=um1(i,j,k,1)*(c2*(um1(i,j,k,5)/um1(i,j,k,1)        &
                -c1*(um1(i,j,k,2)**2+um1(i,j,k,3)**2 &
                +um1(i,j,k,4)**2)/um1(i,j,k,1)**2.))
             pe=max(0.05,pe)
             pres(i,j,k)=pe
             um1(i,j,k,5)=(1.0/c2)*pe+c1*                &
             (um1(i,j,k,2)**2+um1(i,j,k,3)**2+um1(i,j,k,4)**2)/um1(i,j,k,1)
             re=um1(i,j,k,1)
             re=max(0.05,re)
             um1(i,j,k,1)=re
            ENDDO
           ENDDO
          ENDDO
         CALL frontera_x_val (1)
         CALL frontera_z_val (1)
         CALL frontera_y_val (1)
         CALL varc_to_varnc()


         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

        write(6,*)temp(10,1,1),temp(10,1,109)

!______________________________________________________________________
!  2
         write(6,*)'2'
         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         CALL esponja()
         CALL frontera_x()
         CALL frontera_z()
         CALL frontera_y()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx

             um1(i,j,k,m)=0.75*um0(i,j,k,m)+  &
                          0.25*(um1(i,j,k,m)+ &
                          dt*(rs(i,j,k)+rsv(i,j,k)))
            ENDDO
           ENDDO
          ENDDO

          CALL FXSOURCE (2,m)
          CALL FZSOURCE (2,m)
          CALL FYSOURCE (2,m)
         ENDDO
         
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             pe=um1(i,j,k,1)*(c2*(um1(i,j,k,5)/um1(i,j,k,1)        &
                -c1*(um1(i,j,k,2)**2+um1(i,j,k,3)**2 &
                +um1(i,j,k,4)**2)/um1(i,j,k,1)**2.))
             pe=max(0.05,pe)
             pres(i,j,k)=pe
             um1(i,j,k,5)=(1.0/c2)*pe+c1*                &
             (um1(i,j,k,2)**2+um1(i,j,k,3)**2+um1(i,j,k,4)**2)/um1(i,j,k,1)
             re=um1(i,j,k,1)
             re=max(0.05,re)
             um1(i,j,k,1)=re
            ENDDO
           ENDDO
          ENDDO

         CALL frontera_x_val (2)
         CALL frontera_z_val (2)
         CALL frontera_y_val (2)
         CALL varc_to_varnc()
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo
        write(6,*)temp(10,1,1),temp(10,1,109)

!______________________________________________________________________
!  3
         write(6,*)'3'
         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
         CALL esponja()
         CALL frontera_x()
         CALL frontera_z()
         CALL frontera_y()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx

             um1(i,j,k,m)=(1./3.)*(um0(i,j,k,m)+   &
                                   2.*um1(i,j,k,m)+  &
                                   2.*dt*(rs(i,j,k)+rsv(i,j,k)))

            ENDDO
           ENDDO
          ENDDO
          CALL FXSOURCE (3,m)
          CALL FZSOURCE (3,m)
          CALL FYSOURCE (3,m)
         ENDDO

        write(6,*)temp(10,1,1),temp(10,1,109)

         Write(6,*) 'FILTRADO' , itf
         CALL ESTADISTICAS()
!         DO m=1,nd
!         CALL SOURCE(m)
!         ENDDO
         CALL frontera_x_val (3)
         CALL frontera_z_val (3)
         CALL frontera_y_val (3)

         CALL FILTRADO(itf)
         itf=itf+1
         if(itf.eq.7) itf=1

          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             pe=um1(i,j,k,1)*(c2*(um1(i,j,k,5)/um1(i,j,k,1)        &
                -c1*(um1(i,j,k,2)**2+um1(i,j,k,3)**2 &
                +um1(i,j,k,4)**2)/um1(i,j,k,1)**2.))
             pe=max(0.05,pe)
             pres(i,j,k)=pe
             um1(i,j,k,5)=(1.0/c2)*pe+c1*                &
             (um1(i,j,k,2)**2+um1(i,j,k,3)**2+um1(i,j,k,4)**2)/um1(i,j,k,1)
             re=um1(i,j,k,1)
             re=max(0.05,re)
             um1(i,j,k,1)=re
            ENDDO
           ENDDO
          ENDDO

         CALL varc_to_varnc()


         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um0(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

!         write(63,123)dto,u(30,80,54,1),pres(30,80,54),um1(30,80,54,1)
!         write(64,123)dto,u(4,55,71,1),pres(4,55,71),um1(4,55,71,1)
!         write(65,123)dto,u(4,55,77,1),pres(4,55,77),um1(4,55,77,1)
!         write(66,123)dto,u(4,55,81,1),pres(4,55,81),um1(4,55,81,1)
!         write(67,123)dto,u(7,55,81,1),pres(7,55,81),um1(7,55,81,1)
!         write(68,123)dto,u(7,55,77,1),pres(7,55,77),um1(7,55,77,1)
!         write(69,123)dto,u(7,55,73,1),pres(7,55,73),um1(7,55,73,1)
!         write(70,123)dto,u(7,55,70,1),pres(7,55,70),um1(7,55,70,1)
!         write(71,123)dto,u(9,55,71,1),pres(9,55,71),um1(9,55,71,1)
!         write(72,123)dto,u(9,55,75,1),pres(9,55,75),um1(9,55,75,1)
!         write(73,123)dto,u(9,55,79,1),pres(9,55,79),um1(9,55,79,1)
!         write(74,123)dto,u(9,55,81,1),pres(9,55,81),um1(9,55,81,1)
!         write(75,123)dto,u(9,55,73,1),pres(9,55,73),um1(9,55,73,1)
!         write(76,123)dto,u(4,55,75,1),pres(4,55,75),um1(4,55,75,1)
!         write(77,123)dto,u(2,50,69,1),pres(2,50,69),um1(2,50,69,1)
!         write(78,123)dto,u(4,50,71,1),pres(4,50,71),um1(4,50,71,1)
!         write(79,123)dto,u(4,60,77,1),pres(4,60,77),um1(4,60,77,1)
!         write(80,123)dto,u(4,60,81,1),pres(4,60,81),um1(4,60,81,1)
!         write(81,123)dto,u(7,58,77,1),pres(7,58,77),um1(7,58,77,1)
!         write(82,123)dto,u(7,48,73,1),pres(7,48,73),um1(7,48,73,1)
!         write(83,123)dto,u(7,48,70,1),pres(7,48,70),um1(7,48,70,1)
!         write(85,123)dto,u(9,62,71,1),pres(9,62,71),um1(9,62,71,1)
!         write(85,123)dto,u(9,57,75,1),pres(9,57,75),um1(9,57,75,1)
!         write(86,123)dto,u(9,53,79,1),pres(9,53,79),um1(9,53,79,1)
!         write(87,123)dto,u(9,53,81,1),pres(9,53,81),um1(9,53,81,1)
!         write(88,123)dto,u(9,53,73,1),pres(9,53,73),um1(9,53,73,1)

         
        avieu=int(dto/tvieu)
        IF(avieu.ge.1) THEN
         CALL MAXVAL()
         ivieu=ivieu+1
         tvieu=float(ivieu)*tvieu1
        endif

        agrava=int(dto/dgrava)
        IF(agrava.ge.1) THEN
         CALL GRAVA ()
        ENDIF
         
!        write(6,*)temp(10,1,1),temp(10,1,109)
         dto=dto+dt
         it=it+1
        ENDDO
        RETURN
 101    FORMAT(a3,i4,a1,a3,e12.4,a1,a7,e12.4)
 123    FORMAT(f16.8,f16.8,f16.8,f16.8)
        END SUBROUTINE RUNGEK

