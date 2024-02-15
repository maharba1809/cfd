      SUBROUTINE DERIV_VEL()
      use dimensiones
      use derivtools
      use velocidades
      use derivvel
      use consderper
      use consdernper
      use deltas
      IMPLICIT NONE
      integer i,j,k,l,m



!     VELOCIDADES

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=u(i,j,k,1)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dvel(i,j,k,1)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=u(i,j,k,2)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dvel(i,j,k,2)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=u(i,j,k,3)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dvel(i,j,k,3)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=u(i,j,k,1)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         dvel(i,j,k,4)=dv(j)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=u(i,j,k,2)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         dvel(i,j,k,5)=dv(j)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=u(i,j,k,3)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         dvel(i,j,k,6)=dv(j)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=u(i,j,k,1)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dvel(i,j,k,7)=dw(k)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=u(i,j,k,2)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dvel(i,j,k,8)=dw(k)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=u(i,j,k,3)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dvel(i,j,k,9)=dw(k)
        ENDDO
       END DO
      END DO

!    TEMPERATURA

     DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=temp(i,j,k)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dtemp(i,j,k,1)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=temp(i,j,k)
        ENDDO
        call derivper(ny,dv,ayp,byp,cyp,deltay)
        DO j=1,ny
         dtemp(i,j,k,2)=dv(j)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=temp(i,j,k)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dtemp(i,j,k,3)=dw(k)
        ENDDO
       END DO
      END DO

!     CONCENTRACION
     IF(nd.gt.5) THEN
     DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=conc(i,j,k)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         dconc(i,j,k,1)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=conc(i,j,k)
        ENDDO
        call derivper(ny,dv,ayp,byp,cyp,deltay)
        DO j=1,ny
         dconc(i,j,k,2)=dv(j)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=conc(i,j,k)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         dconc(i,j,k,3)=dw(k)
        ENDDO
       END DO
      END DO
      ENDIF

      return
      end subroutine deriv_vel
