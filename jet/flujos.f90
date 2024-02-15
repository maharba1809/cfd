      SUBROUTINE flujos(neq)


      use derivvel
      use jacobtools
      use mflujos
      use viscosidades
      use consadim
      use dimensiones
      use variables
      use velocidades
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: p0,p1,p2
      integer :: neq,mi


      IF (neq.eq.1) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)
            f(i,j,k)=-um(i,j,k,3)
            g(i,j,k)=-um(i,j,k,4)
         END DO
        END DO
       END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            ev(i,j,k)=0.0
            fv(i,j,k)=0.0
            gv(i,j,k)=0.0
         END DO
        END DO
       END DO

      ELSEIF (neq.eq.2) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,1)-c3*pres(i,j,k)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,1)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,1)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           ev(i,j,k)=vis(i,j,k)*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           fv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           gv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
            gv(i,j,1)=0.0
            gv(i,j,nz)=0.0
          END DO
         END DO


      ELSEIF (neq.eq.3) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,2)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,2)-c3*pres(i,j,k)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           ev(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           fv(i,j,k)=vis(i,j,k)*(2.*p1-p0-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           gv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
            gv(i,j,1)=0.0
            gv(i,j,nz)=0.0
          END DO
         END DO


      ELSEIF (neq.eq.4) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           e(i,j,k)=-um(i,j,k,2)*u(i,j,k,3)
           f(i,j,k)=-um(i,j,k,3)*u(i,j,k,3)
           g(i,j,k)=-um(i,j,k,4)*u(i,j,k,3)-c3*pres(i,j,k)
          END DO
         END DO
        END DO
        DO k=2,nz-1
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           ev(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
           ev(i,j,1)=0.0
           ev(i,j,nz)=0.0
          END DO
         END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           fv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
           fv(i,j,1)=0.0
           fv(i,j,nz)=0.0
          END DO
         END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           gv(i,j,k)=vis(i,j,k)*(2.*p2-p1-p0)
          END DO
         END DO
        END DO

      ELSEIF (neq.eq.5) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-u(i,j,k,1)*(um(i,j,k,5)+pres(i,j,k))
            f(i,j,k)=-u(i,j,k,2)*(um(i,j,k,5)+pres(i,j,k))
            g(i,j,k)=-u(i,j,k,3)*(um(i,j,k,5)+pres(i,j,k))
          END DO
         END DO
        END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ev(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,1)*dtemp(i,j,k,1)+    &
                     jbn(i,j,k,4)*dtemp(i,j,k,2)+    &
                     jbn(i,j,k,7)*dtemp(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fv(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,2)*dtemp(i,j,k,1)+    &
                     jbn(i,j,k,5)*dtemp(i,j,k,2)+    &
                     jbn(i,j,k,8)*dtemp(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=2,nz-1
          DO j=1,ny
           DO i=1,nx
            gv(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,3)*dtemp(i,j,k,1)+    &
                     jbn(i,j,k,6)*dtemp(i,j,k,2)+    &
                     jbn(i,j,k,9)*dtemp(i,j,k,3))
           END DO
          END DO
         END DO
          DO j=1,ny
           DO i=1,nx
            gv(i,j,1)=0.0
            gv(i,j,nz)=0.0
           ENDDO
          ENDDO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           ev(i,j,k)=ev(i,j,k)+u(i,j,k,1)*vis2(i,j,k)*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           ev(i,j,k)=ev(i,j,k)+u(i,j,k,2)*vis2(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=2,nz-1
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,3)*dvel(i,j,k,3)+ &
              jbn(i,j,k,6)*dvel(i,j,k,6)+ &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+ &
              jbn(i,j,k,4)*dvel(i,j,k,4)+ &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           ev(i,j,k)=ev(i,j,k)+u(i,j,k,3)*vis2(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           fv(i,j,k)=u(i,j,k,1)*vis2(i,j,k)*(p1+p2)+fv(i,j,k)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           fv(i,j,k)=vis2(i,j,k)*(2.e0*p1-p0-p2)  &
                   *u(i,j,k,2)+fv(i,j,k)
          END DO
         END DO
        END DO
        DO k=2,nz-1
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           fv(i,j,k)=vis2(i,j,k)*(p1+p2)    &
                   *u(i,j,k,3)+fv(i,j,k)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           p2=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           gv(i,j,k)=vis2(i,j,k)*(p1+p2)         &
                   *u(i,j,k,1)+gv(i,j,k)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
           gv(i,j,1)=0.0
           gv(i,j,nz)=0.0
          END DO
         END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           gv(i,j,k)=vis2(i,j,k)*(p1+p2)    &
                   *u(i,j,k,2)+gv(i,j,k)
          END DO
         END DO
        END DO
         DO j=1,ny
          DO i=1,nx
           gv(i,j,1)=0.0
           gv(i,j,nz)=0.0
          END DO
         END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)+  &
              jbn(i,j,k,4)*dvel(i,j,k,4)+  &
              jbn(i,j,k,7)*dvel(i,j,k,7)
           p1=jbn(i,j,k,2)*dvel(i,j,k,2)+  &
              jbn(i,j,k,5)*dvel(i,j,k,5)+  &
              jbn(i,j,k,8)*dvel(i,j,k,8)
           p2=jbn(i,j,k,3)*dvel(i,j,k,3)+  &
              jbn(i,j,k,6)*dvel(i,j,k,6)+  &
              jbn(i,j,k,9)*dvel(i,j,k,9)
           gv(i,j,k)=vis2(i,j,k)*(2.e0*p2-p1-p0) &
                   *u(i,j,k,3)+gv(i,j,k)
          END DO
         END DO
        END DO

      ELSEIF (neq.eq.6) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-u(i,j,k,1)*um(i,j,k,6)
            f(i,j,k)=-u(i,j,k,2)*um(i,j,k,6)
            g(i,j,k)=-u(i,j,k,3)*um(i,j,k,6)
          END DO
         END DO
        END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ev(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,1)*dconc(i,j,k,1)+    &
                     jbn(i,j,k,4)*dconc(i,j,k,2)+    &
                     jbn(i,j,k,7)*dconc(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fv(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,2)*dconc(i,j,k,1)+    &
                     jbn(i,j,k,5)*dconc(i,j,k,2)+    &
                     jbn(i,j,k,8)*dconc(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gv(i,j,k)=vis6(i,j,k)*(          &
                     jbn(i,j,k,3)*dconc(i,j,k,1)+    &
                     jbn(i,j,k,6)*dconc(i,j,k,2)+    &
                     jbn(i,j,k,9)*dconc(i,j,k,3))
           END DO
          END DO
         END DO
      ENDIF
      return
      end subroutine flujos
