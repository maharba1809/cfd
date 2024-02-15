      SUBROUTINE divergencia(neq)
      use dimensiones
      use derivtools
      use mflujos
      use jacobtools
      use dmflujos
      use right
      use consderper
      use consdernper
      use deltas
      use bob
      use variables

      IMPLICIT NONE
      INTEGER :: neq
      integer i,j,k,l,m



         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            e(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,1)*e(i,j,k)+ &
                                    jbn(i,j,k,2)*f(i,j,k)+ &
                                    jbn(i,j,k,3)*g(i,j,k))
            f(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,4)*e(i,j,k)+ &
                                    jbn(i,j,k,5)*f(i,j,k)+ &
                                    jbn(i,j,k,6)*g(i,j,k))
            g(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,7)*e(i,j,k)+ &
                                    jbn(i,j,k,8)*f(i,j,k)+ &
                                    jbn(i,j,k,9)*g(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ev(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,1)*ev(i,j,k)+ &
                                     jbn(i,j,k,2)*fv(i,j,k)+ &
                                     jbn(i,j,k,3)*gv(i,j,k))
            fv(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,4)*ev(i,j,k)+ &
                                     jbn(i,j,k,5)*fv(i,j,k)+ &
                                     jbn(i,j,k,6)*gv(i,j,k))
            gv(i,j,k)=jbn(i,j,k,10)*(jbn(i,j,k,7)*ev(i,j,k)+ &
                                     jbn(i,j,k,8)*fv(i,j,k)+ &
                                     jbn(i,j,k,9)*gv(i,j,k))
           END DO
          END DO
         END DO


      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=e(i,j,k)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dere(i,j,k)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=ev(i,j,k)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         derev(i,j,k)=du(i)
        ENDDO
       END DO
      END DO



      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=f(i,j,k)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         derf(i,j,k)=dv(j)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=fv(i,j,k)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         derfv(i,j,k)=dv(j)
        ENDDO
       END DO
      END DO



      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=g(i,j,k)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         derg(i,j,k)=dw(k)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=gv(i,j,k)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dergv(i,j,k)=dw(k)
        ENDDO
       END DO
      END DO

!FRONTERA CERO FLUJO CONVECTIVO en x=1,nx

      DO k=1,nz
       DO j=1,ny
        dere(1,j,k)=0.0
        dere(nx,j,k)=0.0
        derev(nx,j,k)=0.0
       ENDDO
      ENDDO
      DO j=1,ny
       DO i=1,nx
        derg(i,j,1)=0.0
        derg(i,j,nz)=0.0
       ENDDO
      ENDDO
      DO k=1,nz
       DO i=1,nx
        derf(i,1,k)=0.0
        derf(i,ny,k)=0.0
       ENDDO
      ENDDO
!
!      do i=1,nx
!       write(6,*)i,neq,dere(i,50,50), &
!                       derev(i,50,50)
!       write(6,*)' ',derf(i,50,50)  &
!                    ,derfv(i,50,50)
!       write(6,*)' ',derg(i,50,50)  &
!                    ,dergv(i,50,50)
!       write(6,*)' '
!      enddo

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         rs(i,j,k)=jbn(i,j,k,11)*(dere(i,j,k)+derf(i,j,k)+derg(i,j,k)) &
                   -esource(i)*(um(i,j,k,neq)-esplayer(j,k,neq))
         rsv(i,j,k)=jbn(i,j,k,11)*(derev(i,j,k)+derfv(i,j,k)+dergv(i,j,k))
        END DO
       END DO
      END DO



      end subroutine divergencia
