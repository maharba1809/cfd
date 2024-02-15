      SUBROUTINE JACOBEANO()
      use dimensiones
      use jacobtools
      use derivtools
      use mallagrid
      use consderper
      use consdernper
      use deltas
      IMPLICIT NONE
      integer :: i,j,k,l,m



      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=x(i,j,k,1)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         jbn(i,j,k,1)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=x(i,j,k,2)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         jbn(i,j,k,4)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         du(i)=x(i,j,k,3)
        ENDDO
        call derivnper(nx,du,axnp,bxnp,cxnp,deltax)
        DO i=1,nx
         jbn(i,j,k,7)=du(i)
        ENDDO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=x(i,j,k,1)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         jbn(i,j,k,2)=dv(j)
        ENDDO
       END DO
      END DO


      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=x(i,j,k,2)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         jbn(i,j,k,5)=dv(j)
        ENDDO
       END DO
      END DO


      DO k=1,nz
       DO i=1,nx
        DO j=1,ny
         dv(j)=x(i,j,k,3)
        ENDDO
        call derivnper(ny,dv,aynp,bynp,cynp,deltay)
        DO j=1,ny
         jbn(i,j,k,8)=dv(j)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=x(i,j,k,1)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         jbn(i,j,k,3)=dw(k)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=x(i,j,k,2)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         jbn(i,j,k,6)=dw(k)
        ENDDO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        DO k=1,nz
         dw(k)=x(i,j,k,3)
        ENDDO
        call derivnper(nz,dw,aznp,bznp,cznp,deltaz)
        DO k=1,nz
         jbn(i,j,k,9)=dw(k)
        ENDDO
       END DO
      END DO


!      Inversion de matriz

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbn(i,j,k,10)=jbn(i,j,k,1)*jbn(i,j,k,5)*jbn(i,j,k,9) &
                    +jbn(i,j,k,2)*jbn(i,j,k,6)*jbn(i,j,k,7)  & 
                    +jbn(i,j,k,3)*jbn(i,j,k,4)*jbn(i,j,k,8)  &
                    -jbn(i,j,k,1)*jbn(i,j,k,6)*jbn(i,j,k,8)  &
                    -jbn(i,j,k,2)*jbn(i,j,k,4)*jbn(i,j,k,9)  &
                    -jbn(i,j,k,3)*jbn(i,j,k,5)*jbn(i,j,k,7)
         If(abs(jbn(i,j,k,10)).le.1.e-20) THEN
          print *,'  VALOR IGUAL A CERO EN :', i,j,k
          STOP
         ENDIF

        END DO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbn(i,j,k,11)=1./jbn(i,j,k,10)
         t1=jbn(i,j,k,5)*jbn(i,j,k,9)-jbn(i,j,k,6)*jbn(i,j,k,8)
         t2=jbn(i,j,k,4)*jbn(i,j,k,9)-jbn(i,j,k,6)*jbn(i,j,k,7)
         t3=jbn(i,j,k,4)*jbn(i,j,k,8)-jbn(i,j,k,5)*jbn(i,j,k,7)
         t4=jbn(i,j,k,2)*jbn(i,j,k,9)-jbn(i,j,k,3)*jbn(i,j,k,8)
         t5=jbn(i,j,k,1)*jbn(i,j,k,9)-jbn(i,j,k,3)*jbn(i,j,k,7)
         t6=jbn(i,j,k,1)*jbn(i,j,k,8)-jbn(i,j,k,2)*jbn(i,j,k,7)
         t7=jbn(i,j,k,2)*jbn(i,j,k,6)-jbn(i,j,k,3)*jbn(i,j,k,5)
         t8=jbn(i,j,k,1)*jbn(i,j,k,6)-jbn(i,j,k,3)*jbn(i,j,k,4)
         t9=jbn(i,j,k,1)*jbn(i,j,k,5)-jbn(i,j,k,2)*jbn(i,j,k,4)
         jbn(i,j,k,1)= t1*jbn(i,j,k,11)
         jbn(i,j,k,2)=-t2*jbn(i,j,k,11)
         jbn(i,j,k,3)= t3*jbn(i,j,k,11)
         jbn(i,j,k,4)=-t4*jbn(i,j,k,11)
         jbn(i,j,k,5)= t5*jbn(i,j,k,11)
         jbn(i,j,k,6)=-t6*jbn(i,j,k,11)
         jbn(i,j,k,7)= t7*jbn(i,j,k,11)
         jbn(i,j,k,8)=-t8*jbn(i,j,k,11)
         jbn(i,j,k,9)= t9*jbn(i,j,k,11)
        END DO
       END DO
      END DO

      end subroutine jacobeano
