!-----------------------------------------------------------------------
      program compacto
!-----------------------------------------------------------------------

      use dimensiones
    

      nx=100
      ny=100
      nz=100
      call flow_alloc()
      call inideriv()
      call deriva()

      end program compacto


      SUBROUTINE deriva()
      use dimensiones 
      use consderper
      use consdernper
      real, dimension (nx) ::  u
      real :: pi,delta


      pi=2.*acos(0.)
      delta=2.*Pi/float(nx)
      do i=1,nx
       u(i)=sin(float(i)*delta)
      enddo
      delta=1./delta
      
      call derivper(nx,u,axp,bxp,cxp,delta)
!      call derivnper(nx,u,axnp,bxnp,cxnp,delta)

      return
      end subroutine deriva
