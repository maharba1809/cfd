
        SUBROUTINE inideriv()
        use cderivper
        use cderivnper
        use mderivper
        use indices
  
!       NO PERIODICO

        DO i=1,nx
         bnpx(i)=1.
        ENDDO
        DO i=1,ny
         bnpy(i)=1.
        ENDDO
        DO i=1,nz
         bnpz(i)=1.
        ENDDO

        anpx(1)=0.0
        anpx(2)=0.25
        anpx(nx-1)=0.25
        anpx(nx)=2.
        DO i=3,nx-2
         apx(i)=1./3.
        ENDDO
        anpy(1)=0.0
        anpy(2)=0.25
        anpy(ny-1)=0.25
        anpy(n)=2.
        DO i=3,ny-2
         anpy(i)=1./3.
        ENDDO
        anpz(1)=0.0
        anpz(2)=0.25
        anpz(nz-1)=0.25
        anpz(nz)=2.
        DO i=3,nz-2
         anpz(i)=1./3.
        ENDDO

        
        cnpx(1)=2.0
        cnpx(2)=0.25
        cnpx(nx-1)=0.25
        cnpx(nx)=0.0
        DO i=3,nx-2
         cnpx(i)=1./3.
        ENDDO
        cnpy(1)=2.0
        cnpy(2)=0.25
        cnpy(ny-1)=0.25
        cnpy(ny)=0.0
        DO i=3,ny-2
         cnpy(i)=1./3.
        ENDDO
        cnpz(1)=2.0
        cnpz(2)=0.25
        cnpz(nz-1)=0.25
        cnpz(nz)=0.0
        DO i=3,nz-1
         cnpz(i)=1./3.
        ENDDO
        
        ars=7./9.
        brs=1./36.
        ars1=-5./2.
        brs1=2.
        crs1=1./2.
        ars2=3./4.


!        PERIODICO
        DO i=1,nx
         bpx(i)=1.
         apx(i)=1./3.
         cpx(i)=1./3.
        ENDDO
        DO i=1,ny
         bpy(i)=1.
         apy(i)=1./3.
         cpy(i)=1./3.
        ENDDO
        DO i=1,nz
         bpz(i)=1.
         apz(i)=1./3.
         cpz(i)=1./3.
        ENDDO


        SUBROUTINE derivnper(n,u,ars,brs,crs,delta)

        IMPLICIT NONE
        real ars(n),brs(n),crs(n)
        real bet
        use cderivnper
        use mderivnper
        use indices
        

        DO i=3,n-2
         rs(i)=delta*(als*(u(i+1)-u(i-1))+
     1               bls*(u(i+2)-u(i-2)))
        ENDDO
         rs(1)=delta*(als1*u(1)+bls1*u(2)+cls1(3))
         rs(n)=delta*(-als1*u(n)-bls1*u(n-1)-cls1(n-2))
         rs(2)=delta*(als2*(u(3)-u(1))
         rs(n-1)=delta*(als2*(u(n)-u(n-2))

        bet=brs(1)
        u(1)=r(1)/bet
!
        do j=2,n
         gam(j)=crs(j-1)/bet
         bet=brs(j)-ars(j)*gam(j)
         u(j)=(r(j)-ars(j)*u(j-1))/bet
        enddo
!
        do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
        enddo
        RETURN
        END

        SUBROUTINE derivper(n,u,ars,brs,crs,delta)

        IMPLICIT NONE
        real delta,bet
        real ars(n),brs(n),crs(n)
        integer n
        use cderivper
        use mderivper
        use indices


        DO i=3,n-2
         rs(i)=delta*(als*(u(i+1)-u(i-1))+    &
                    bls*(u(i+2)-u(i-2)))
        ENDDO
         rs(1)=delta*(als*(u(2)-u(n))+    &
                      bls*(u(3)-u(n-1)))
         rs(n)=delta*(als*(u(1)-u(n-1))+    &
                    bls*(u(2)-u(n-2)))
         rs(2)=delta*(als*(u(3)-u(1))+    &
                    bls*(u(4)-u(n)))
         rs(n-1)=delta*(als*(u(n)-u(n-2))+    &
                    bls*(u(i+2)-u(i-2)))

        alfa=a(1)
        beta=b(1)
        gamma=-b(1)
        bb(1)=b(1)-gamma
        bb(n)=b(n)-alfa*beta/gamma
        do i=1,n-1
         bb(i)=b(i)
        enddo
!       PRIMERA TRIDIAGONAL

        bet=bb(1)
        x(1)=r(1)/bet

        do j=2,n
         gam(j)=crs(j-1)/bet
         bet=bb(j)-ars(j)*gam(j)
         x(j)=(r(j)-ars(j)*x(j-1))/bet
        enddo

        do j=n-1,1,-1
         x(j)=x(j)-gam(j+1)*x(j+1)
        enddo
!       ******

        u(1)=gamma
        u(n)=alpha
        do i=2,n-1
         u(i)=0.0
        enddo

!       SEGUNDA TRIDIAGONAL

        bet=bb(1)
        z(1)=u(1)/bet

        do j=2,n
         gam(j)=crs(j-1)/bet
         bet=bb(j)-ars(j)*gam(j)
         z(j)=(u(j)-ars(j)*z(j-1))/bet
        enddo

        do j=n-1,1,-1
         z(j)=z(j)-gam(j+1)*z(j+1)
        enddo
!       ******

        fact=(x(1)+beta*x(n)/gamma)/(1+z(1)+beta*z(n)/gamma)
        do i=1,n
         u(i)=x(i)-fact*z(i)
        enddo

        RETURN
        END



