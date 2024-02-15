      use dimensiones
      use consderper
      use consdernper
      use consderper_f
      use consdernper_f
      use consrs
      use consrs_f
      use funciones
      Implicit NONE
      real pi,delta,deltax,deltay,deltaz
      integer i


      CALL INI ()
      CALL MEMORIA ()
      write(6,*)'aqui'
      CALL INIDERIV ()
      CALL INIFILTRO ()
      pi=2.*acos(0.)
      delta=2*Pi/float(nx)
      deltax=1./delta
      deltay=deltax
      deltaz=deltax
      do i=1,nx
      u(i)=sin(float(i)*delta)+sin(30.*float(i)*delta)
      u2(i)=sin(float(i)*delta)+sin(30.*float(i)*delta)
      write(16,100)i,u(i),u2(i)
      enddo
  
!      CALL derivnper(nz,u,aznp,bznp,cznp,deltaz)
!      CALL derivper(nz,u2,azp,bzp,czp,deltaz)
      CALL filtronper(nz,u,aznpf,bznpf,cznpf)
      CALL filtroper(nz,u2,azpf,bzpf,czpf)

      write(6,*)'DERIVADA'
      do i=1,nx
      write(17,100)i,u(i),u2(i)
      enddo
      CALL MEMORIA_DEALLOC ()
 100  FORMAT(i3,e16.6,e16.6)

      end

      SUBROUTINE ini()

      use dimensiones
      use consfiltro

      nx=100.
      ny=100.
      nz=100.
      cfilt=0.49

      return
      end subroutine ini
