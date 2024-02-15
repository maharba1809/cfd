!________________________________________________________
      SUBROUTINE INIESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      IMPLICIT NONE
      integer i,j,k,l,m

       do k=1,nz
       do j=1,ny
       do i=1,nx
       do l=1,15
        st(i,j,k,l)=0.0
       enddo
       enddo
       enddo
       enddo
       RETURN
       END SUBROUTINE INIESTADISTICAS


!________________________________________________________
      SUBROUTINE ESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      use tiempo
      IMPLICIT NONE

      integer  i,j,k,l,m


      do k=1,nz
      do j=1,ny
      do i=1,nx
       st(i,j,k,1)=st(i,j,k,1)+dt*u(i,j,k,1)
       st(i,j,k,2)=st(i,j,k,2)+dt*u(i,j,k,1)**2.
       st(i,j,k,3)=st(i,j,k,3)+dt*u(i,j,k,2)
       st(i,j,k,4)=st(i,j,k,4)+dt*u(i,j,k,2)**2.
       st(i,j,k,5)=st(i,j,k,5)+dt*u(i,j,k,3)
       st(i,j,k,6)=st(i,j,k,6)+dt*(u(i,j,k,3)*u(i,j,k,3))
       st(i,j,k,7)=st(i,j,k,7)+dt*(u(i,j,k,1)*u(i,j,k,2))
       st(i,j,k,8)=st(i,j,k,8)+dt*(u(i,j,k,1)*u(i,j,k,3))
       st(i,j,k,9)=st(i,j,k,9)+dt*(u(i,j,k,2)*u(i,j,k,3))
       st(i,j,k,10)=st(i,j,k,10)+dt*pres(i,j,k)
       st(i,j,k,11)=st(i,j,k,11)+dt*pres(i,j,k)**2
       st(i,j,k,12)=st(i,j,k,12)+dt*um1(i,j,k,1)
       st(i,j,k,13)=st(i,j,k,13)+dt*um1(i,j,k,1)**2
       st(i,j,k,14)=st(i,j,k,14)+dt*temp(i,j,k)
       st(i,j,k,15)=st(i,j,k,15)+dt*temp(i,j,k)**2
      enddo
      enddo
      enddo

     ! write(57,100)dto,um1(10,55,55,1),um1(20,55,55,1),um1(30,55,55,1),um1(40,55,55,1),um1(50,55,55,1)
!      write(58,100)dto,pres(10,55,55),pres(20,55,55),pres(30,55,55),pres(40,55,55),pres(50,55,55)
!      write(59,100)dto,pres(10,40,40),pres(20,40,40),pres(30,40,40),pres(40,40,40),pres(50,40,40)
      !do j=1,ny
      !   write(*,*)pres(1,j,nz/2),u(1,j,nz/2,1)
      !enddo
      
      write(60,100)dto,pres(5,ny/2+10,nz/2+10)
      write(61,100)dto,pres(5,ny/2+20,nz/2+20)
      write(62,100)dto,pres(10,ny/2+10,nz/2+10)
      write(63,100)dto,pres(10,ny/2+20,nz/2+20)
      write(64,100)dto,pres(10,ny/2+30,nz/2+30)
      write(65,100)dto,pres(20,ny/2+10,nz/2+10)
      write(66,100)dto,pres(20,ny/2+20,nz/2+20)
      write(67,100)dto,pres(20,ny/2+30,nz/2+30)
      write(68,100)dto,pres(40,ny/2+20,nz/2+20)
      write(69,100)dto,pres(40,ny/2+30,nz/2+30)
      RETURN
 100  FORMAT(f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)
      END SUBROUTINE ESTADISTICAS

!________________________________________________________
      SUBROUTINE FINESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      use tiempo
      IMPLICIT NONE
      integer i,j,k,l,m

       do k=1,nz
       do j=1,ny
       do i=1,nx
       do l=1,15
        st(i,j,k,l)=st(i,j,k,l)/(dto)
       enddo
       enddo
       enddo
       enddo

       !write(6,*)st(100,55,55,1),st(100,55,55,2)
       OPEN(67,file='stat.data',form='unformatted')
       write(67)nx,ny,nz
       write(67)st
       CLOSE(67)
       RETURN
       END SUBROUTINE FINESTADISTICAS


