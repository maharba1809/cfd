module vary
integer i,j,k,nx,ny,nz,cont,ini,fin,paso,di,in_Nvar,sample
real lx,ly,lz,rho,ma
real,allocatable,dimension(:,:,:)::u,v,w,p,t,x,y,z
character(3)::campo
character(2)::serie
end module
!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
program tpmaker
use vary
implicit none
!call leer_datos(nx,ny,nz,lx,ly,lz)


open(66,file='grid.bin',form='unformatted')
    read(66)in_Nvar
    read(66)nx,ny,nz
    
    call aloca()
    
    read(66)x
    read(66)y
    read(66)z
close(66)

write(*,*)nx,ny,nz
write(*,*)'initial file:' 
read(*,*)ini  
write(*,*)'final file:' 
read(*,*)fin 
!ini = 0
!fin = 8
paso=1
write(*,*)ini,fin,paso

write(*,*)'serie:' 
!serie='04'
read(*,*)serie

do cont=ini,fin,paso
   write(campo,300)cont
300  format(i3.3)
   write(*,*)'field_'//serie//'.'//campo//''
   open(19,file='field_'//serie//'.'//campo//'',form='unformatted')
       READ(19)in_Nvar
       READ(19)nx,ny,nz
       READ(19)u
       READ(19)v
       READ(19)w
       READ(19)t
       READ(19)p
   close(19)
   open(234,file='jet_'//serie//'_'//campo//'.dat')
     write(*,*)'Escribiendo ---> jet_'//serie//'_'//campo//'.dat'
     !write(234,*)'Variables= "X", "Y", "Z", "U", "V", "W", "P","T","Rho","Ma"'
     write(234,*)'Variables= "X", "Y", "U", "V", "P","T","RHO","MACH"'
     !write(234,*)'Zone I= ',nx,' J= ',ny,'K= ',nz,' Datapacking=point'
     write(234,*)'Zone I= ',nx,' J= ',ny,' Datapacking=point'
     
     do k=nz/2,nz/2
        do j=1,ny
           do i=1,nx
            rho=p(i,j,k)/t(i,j,k)
            ma = ((u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)**0.5)/(t(i,j,k)**0.5) *0.2
            !write(234,*)x(i,j,k),y(i,j,k),z(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),t(i,j,k),rho,ma
            write(234,*)x(i,j,k),y(i,j,k),u(i,j,k),v(i,j,k),p(i,j,k),t(i,j,k),rho,ma
            if (i.eq.1)then
               write(*,*)j,p(i,j,k),rho
            end if
           end do
         end do
      end do
   close(234)
end do
call dealoca()

end

!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
subroutine leer_datos(nx,ny,nz,lx,ly,lz)
implicit none
integer nx,ny,nz
real lx,ly,lz

open(32,file='data.in',form='formatted')
read(32,100)nx
read(32,100)ny
read(32,100)nz
read(32,200)lx
read(32,200)ly
read(32,200)lz
close(32)

100  format(16x,i4)
200  format(16x,f10.4)
return
end subroutine
!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
subroutine aloca()
use vary

allocate(x(nx,ny,nz))
allocate(y(nx,ny,nz))
allocate(z(nx,ny,nz))
allocate(u(nx,ny,nz))
allocate(v(nx,ny,nz))
allocate(w(nx,ny,nz))
allocate(p(nx,ny,nz))
allocate(t(nx,ny,nz))

return
end subroutine
!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
subroutine dealoca()
use vary

deallocate(x)
deallocate(y)
deallocate(z)
deallocate(u)
deallocate(v)
deallocate(w)
deallocate(p)
deallocate(t)

return
end subroutine
!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
