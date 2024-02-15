
      SUBROUTINE inifiltro()
      use consderper_f
      use consdernper_f
      use consrs_f
      use dimensiones
      use consfiltro
      IMPLICIT NONE
      integer i,j,k,l,m
      real cfilt2,cfilt3,cfilt4,cfilt5

      do i=1,nx
       bxnpf(i)=1.0
      enddo
      do i=1,ny
       bynpf(i)=1.0
      enddo
      do i=1,nz
       bznpf(i)=1.0
      enddo

      axnpf(1)=0.0
      axnpf(nx)=0.0
      do i=2,nx-2
       axnpf(i)=cfilt
      enddo
      aynpf(1)=0.0
      aynpf(ny)=0.0
      do i=2,ny-1
       aynpf(i)=cfilt
      enddo
      aznpf(1)=0.0
      aznpf(nz)=0.0
      do i=2,nz-1
       aznpf(i)=cfilt
      enddo

      cxnpf(1)=0.0
      cxnpf(nx)=0.0
      do i=2,nx-1
       cxnpf(i)=cfilt
      enddo
      aynpf(1)=0.0
      aynpf(ny)=0.0
      do i=2,ny-1
       cynpf(i)=cfilt
      enddo
      aznpf(1)=0.0
      aznpf(nz)=0.0
      do i=2,nz-1
       cznpf(i)=cfilt
      enddo

      do i=1,nx
       bxpf(i)=1.0
       axpf(i)=cfilt
       cxpf(i)=cfilt
      enddo
      do i=1,ny
       bypf(i)=1.0
       aypf(i)=cfilt
       cypf(i)=cfilt
      enddo
      do i=1,nz
       bzpf(i)=1.0
       azpf(i)=cfilt
       czpf(i)=cfilt
      enddo

     cfilt2=cfilt
     cfilt3=cfilt
     cfilt4=cfilt
     cfilt5=cfilt

     arsf=(193.+126.*cfilt)/256.
     brsf=(105.+302.*cfilt)/256.
     crsf=(-15.+30.*cfilt)/64.
     drsf=(45.-90.*cfilt)/512.
     ersf=(-5.+10.*cfilt)/256.
     frsf=(1.-2*cfilt)/512.

     arsf1=(93.+70.*cfilt2)/128.
     brsf1=(7.+18.*cfilt2)/16.
     crsf1=(-7.+14.*cfilt2)/32.
     drsf1=(1.-2.*cfilt2)/16.
     ersf1=(-1.+2.*cfilt2)/128.

     arsf2=11./16.+(5.*cfilt3)/8.
     brsf2=15./32.+(17.*cfilt3)/16.
     crsf2=-3./16.+(3.*cfilt3)/8.
     drsf2=1./32.-(cfilt3)/16.

     arsf3=5./8.+3.*(cfilt4)/4.
     brsf3=1./2.+cfilt4
     crsf3=-1./8.+(cfilt4)/4.

     arsf4=1./2.+cfilt5
     brsf4=1./2.+cfilt5

     return
     end subroutine inifiltro


      SUBROUTINE filtroper(n,u,a,b,c)
      use consrs_f
      use consfiltro
      IMPLICIT NONE
      integer i,j,k,l,m
      integer :: n
      real :: delta,bet,alpha,beta,gamma,fact
      real,dimension (n) ::  u,r,gam,a,b,c,bb
      real, dimension (n) :: x,z


        DO l=6,n-5
         r(l)=0.5*(arsf*(u(l)+u(l))+     &
                   brsf*(u(l+1)+u(l-1))+ &
                   crsf*(u(l+2)+u(l-2))+ &
                   drsf*(u(l+3)+u(l-3))+ &
                   ersf*(u(l+4)+u(l-4))+ &
                   frsf*(u(l+5)+u(l-5)))
        ENDDO
         r(5)=0.5*(arsf*(u(5)+u(5))+     &
                   brsf*(u(6)+u(4))+     &
                   crsf*(u(7)+u(3))+     &
                   drsf*(u(8)+u(2))+     &
                   ersf*(u(9)+u(1))+   &
                   frsf*(u(10)+u(n)))
         r(4)=0.5*(arsf*(u(4)+u(4))+     &
                   brsf*(u(5)+u(3))+     &
                   crsf*(u(6)+u(2))+     &
                   drsf*(u(7)+u(1))+     &
                   ersf*(u(8)+u(n))+   &
                   frsf*(u(9)+u(n-1)))
         r(3)=0.5*(arsf*(u(3)+u(3))+     &
                   brsf*(u(4)+u(2))+     &
                   crsf*(u(5)+u(1))+     &
                   drsf*(u(6)+u(n))+     &
                   ersf*(u(7)+u(n-1))+   &
                   frsf*(u(8)+u(n-2)))
         r(2)=0.5*(arsf*(u(2)+u(2))+     &
                   brsf*(u(3)+u(1))+     &
                   crsf*(u(4)+u(n))+     &
                   drsf*(u(5)+u(n-1))+   &
                   ersf*(u(6)+u(n-2))+   &
                   frsf*(u(7)+u(n-3)))
         r(1)=0.5*(arsf*(u(1)+u(1))+     &
                   brsf*(u(2)+u(n))+     &
                   crsf*(u(3)+u(n-1))+   &
                   drsf*(u(4)+u(n-2))+   &
                   ersf*(u(5)+u(n-3))+   &
                   frsf*(u(6)+u(n-4)))   
         r(n)=0.5*(arsf*(u(n)+u(n))+     &
                   brsf*(u(1)+u(n-1))+   &
                   crsf*(u(2)+u(n-2))+   &
                   drsf*(u(3)+u(n-3))+   &
                   ersf*(u(4)+u(n-4))+   &
                   frsf*(u(5)+u(n-5)))   
         r(n-1)=0.5*(arsf*(u(n-1)+u(n-1))+ &
                     brsf*(u(n)+u(n-2))+   &
                     crsf*(u(1)+u(n-3))+   &
                     drsf*(u(2)+u(n-4))+   &
                     ersf*(u(3)+u(n-5))+   &
                     frsf*(u(4)+u(n-6)))  
         r(n-2)=0.5*(arsf*(u(n-2)+u(n-2))+ &
                     brsf*(u(n-1)+u(n-3))+ &
                     crsf*(u(n)+u(n-4))+   &
                     drsf*(u(1)+u(n-5))+   &
                     ersf*(u(2)+u(n-6))+   &
                     frsf*(u(3)+u(n-7)))   
         r(n-3)=0.5*(arsf*(u(n-3)+u(n-3))+ &
                     brsf*(u(n-2)+u(n-4))+ &
                     crsf*(u(n-1)+u(n-5))+   &
                     drsf*(u(n)+u(n-6))+   &
                     ersf*(u(1)+u(n-7))+   &
                     frsf*(u(2)+u(n-8)))   
         r(n-4)=0.5*(arsf*(u(n-4)+u(n-4))+ &
                     brsf*(u(n-3)+u(n-5))+ &
                     crsf*(u(n-2)+u(n-6))+   &
                     drsf*(u(n-1)+u(n-7))+   &
                     ersf*(u(n)+u(n-8))+   &
                     frsf*(u(1)+u(n-9)))



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
        end subroutine filtroper



      SUBROUTINE filtronper(n,u,a,b,c)
      use consrs_f
      IMPLICIT NONE
      integer i,j,k,l,m
      integer :: n
      real :: delta,bet
      real, dimension (n) ::  u,r,gam,a,b,c

! O(10)
        DO l=6,n-5
         r(l)=0.5*(arsf*(u(l)+u(l))+     &
                   brsf*(u(l+1)+u(l-1))+ &
                   crsf*(u(l+2)+u(l-2))+ &
                   drsf*(u(l+3)+u(l-3))+ &
                   ersf*(u(l+4)+u(l-4))+ &
                   frsf*(u(l+5)+u(l-5)))
        ENDDO
! O(8)
         r(5)=0.5*(arsf1*(u(5)+u(5))+     &
                   brsf1*(u(5+1)+u(5-1))+ &
                   crsf1*(u(5+2)+u(5-2))+ &
                   drsf1*(u(5+3)+u(5-3))+ &
                   ersf1*(u(5+4)+u(5-4)))
         r(n-4)=0.5*(arsf1*(u((n-4))+u((n-4)))+     &
                   brsf1*(u((n-4)+1)+u((n-4)-1))+   &
                   crsf1*(u((n-4)+2)+u((n-4)-2))+   &
                   drsf1*(u((n-4)+3)+u((n-4)-3))+   &
                   ersf1*(u((n-4)+4)+u((n-4)-4)))
! O(6)

         r(4)=0.5*(arsf2*(u(4)+u(4))+     &
                   brsf2*(u(4+1)+u(4-1))+ &
                   crsf2*(u(4+2)+u(4-2))+ &
                   drsf2*(u(4+3)+u(4-3))) 
         r(n-3)=0.5*(arsf2*(u((n-3))+u((n-3)))+     &
                   brsf2*(u((n-3)+1)+u((n-3)-1))+   &
                   crsf2*(u((n-3)+2)+u((n-3)-2))+   &
                   drsf2*(u((n-3)+3)+u((n-3)-3))) 
!O(4)
         r(3)=0.5*(arsf3*(u(3)+u(3))+     &
                   brsf3*(u(2)+u(4))+     &
                   crsf3*(u(1)+u(5)))     
         r(n-2)=0.5*(arsf3*(u(n-2)+u(n-2))+ &
                     brsf3*(u(n-3)+u(n-1))+ &
                     crsf3*(u(n-4)+u(n)))
!O(2)
         r(2)=0.5*(arsf4*(u(2)+u(2))+       &
                   brsf4*(u(1)+u(3)))     
         r(n-1)=0.5*(arsf4*(u(n-1)+u(n-1))+ &
                     brsf4*(u(n-2)+u(n)))  
! SIN FILTRO

         r(1)=u(1)
         r(n)=u(n)


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
!       do l=1,n
!        write(6,*)l,r(l),u(l)
!       enddo
!       stop


        return
        end subroutine filtronper



