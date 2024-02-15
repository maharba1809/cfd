
      SUBROUTINE inideriv()
      use consderper
      use consdernper
      use consrs
      use dimensiones
      IMPLICIT NONE
      integer i,j,k,l,m

      do i=1,nx
       bxnp(i)=1.0
      enddo
      do i=1,ny
       bynp(i)=1.0
      enddo
      do i=1,nz
       bznp(i)=1.0
      enddo

      axnp(1)=0.0
      axnp(2)=0.25
      axnp(nx-1)=0.25
      axnp(nx)=2.0
      do i=3,nx-2
       axnp(i)=1.0/3.0
      enddo
      aynp(1)=0.0
      aynp(2)=0.25
      aynp(ny-1)=0.25
      aynp(ny)=2.0
      do i=3,ny-2
       aynp(i)=1.0/3.0
      enddo
      aznp(1)=0.0
      aznp(2)=0.25
      aznp(nz-1)=0.25
      aznp(nz)=2.0
      do i=3,nz-2
       aznp(i)=1.0/3.0
      enddo

      cxnp(1)=2.0
      cxnp(2)=0.25
      cxnp(nx-1)=0.25
      cxnp(nx)=0.0
      do i=3,nx-2
       cxnp(i)=1.0/3.0
      enddo
      cynp(1)=2.0
      cynp(2)=0.25
      cynp(ny-1)=0.25
      cynp(ny)=0.0
      do i=3,ny-2
       cynp(i)=1.0/3.0
      enddo
      cznp(1)=2.0
      cznp(2)=0.25
      cznp(nz-1)=0.25
      cznp(nz)=0.0
      do i=3,nz-2
       cznp(i)=1.0/3.0
      enddo

      do i=1,nx
       bxp(i)=1.0
       axp(i)=1.0/3.0
       cxp(i)=1.0/3.0
      enddo
      do i=1,ny
       byp(i)=1.0
       ayp(i)=1.0/3.0
       cyp(i)=1.0/3.0
      enddo
      do i=1,nz
       bzp(i)=1.0
       azp(i)=1.0/3.0
       czp(i)=1.0/3.0
      enddo


     ars=7./9.
     brs=1./36.
     ars1=-2.5
     brs1=2.
     crs1=0.5 
     ars2=0.75


     return
     end subroutine inideriv


      SUBROUTINE derivper(n,u,a,b,c,delta)
      use consrs
      IMPLICIT NONE
      integer :: n
      integer i,j,k,l,m
      real :: delta,bet,alpha,beta,gamma,fact
      real,dimension (n) ::  u,r,gam,a,b,c,bb
      real, dimension (n) :: x,z


        DO l=3,n-2
         r(l)=delta*(ars*(u(l+1)-u(l-1))+    &
                   brs*(u(l+2)-u(l-2)))
        ENDDO
         r(1)=delta*(ars*(u(2)-u(n))+        &
                   brs*(u(3)-u(n-1)))
         r(n)=delta*(ars*(u(1)-u(n-1))+      &
                   brs*(u(2)-u(n-2)))
         r(2)=delta*(ars*(u(3)-u(1))+        &
                   brs*(u(4)-u(n)))
         r(n-1)=delta*(ars*(u(n)-u(n-2))+    &
                   brs*(u(1)-u(n-3)))
        alpha=a(1)
        beta=c(n)
        gamma=-b(1)
        bb(1)=b(1)-gamma
        bb(n)=b(n)-alpha*beta/gamma
        do l=2,n-1
         bb(l)=b(l)
        enddo


        bet=bb(1)
        x(1)=r(1)/bet
        

        do l=2,n
         gam(l)=c(l-1)/bet
         bet=bb(l)-a(l)*gam(l)
         x(l)=(r(l)-a(l)*x(l-1))/bet
        enddo

        do l=n-1,1,-1
         x(l)=x(l)-gam(l+1)*x(l+1)
        enddo

        
        u(1)=gamma
        u(n)=alpha
        do l=2,n-1
         u(l)=0.0
        enddo

        bet=bb(1)
        z(1)=u(1)/bet

        do l=2,n
         gam(l)=c(l-1)/bet
         bet=bb(l)-a(l)*gam(l)
         z(l)=(u(l)-a(l)*z(l-1))/bet
        enddo

        do l=n-1,1,-1
         z(l)=z(l)-gam(l+1)*z(l+1)
        enddo

        fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      
        do l=1,n
         u(l)=x(l)-fact*z(l)
        enddo

        return
        end subroutine derivper


      SUBROUTINE derivnper(n,u,a,b,c,delta)
      use consrs
      IMPLICIT NONE
      integer i,j,k,l,m
      integer :: n
      real :: delta,bet
      real, dimension (n) ::  u,r,gam,a,b,c


        DO l=3,n-2
         r(l)=delta*(ars*(u(l+1)-u(l-1))+     &
                   brs*(u(l+2)-u(l-2)))
        ENDDO
         r(1)=delta*(ars1*u(1)+brs1*u(2)+crs1*u(3))
         r(n)=delta*(-ars1*u(n)-brs1*u(n-1)-crs1*u(n-2))
         r(2)=delta*(ars2*(u(3)-u(1)))
         r(n-1)=delta*(ars2*(u(n)-u(n-2)))

        bet=b(1)
        u(1)=r(1)/bet

        do l=2,n
         gam(l)=c(l-1)/bet
         bet=b(l)-a(l)*gam(l)
         u(l)=(r(l)-a(l)*u(l-1))/bet
        enddo

        do l=n-1,1,-1
         u(l)=u(l)-gam(l+1)*u(l+1)
        enddo

        return
        end subroutine derivnper



